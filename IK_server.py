#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from cv2 import sqrt
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np
import time


def rot_x(q):
	R_x = Matrix([[1, 0, 0],
			[0, cos(q), -sin(q)],
			[0, sin(q), cos(q)]])

	return R_x


def rot_y(q):
	R_y = Matrix([[cos(q), 0, sin(q)],
				[0, 1, 0],
				[-sin(q), 0, cos(q)]])

	return R_y


def rot_z(q):
	R_z = Matrix([[cos(q), -sin(q), 0],
				[sin(q), cos(q), 0],
				[0, 0, 1]])

	return R_z

def forward_kinematics(theta):
	q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
	d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
	a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
	alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

	s = {alpha0:     0.0, a0:      0, d1:  0.75, q1: theta[0],
		alpha1: -pi / 2, a1:   0.35, d2:     0, q2: theta[1] - pi / 2,
		alpha2:     0.0, a2:   1.25, d3:     0, q3: theta[2],
		alpha3: -pi / 2, a3: -0.054, d4:  1.50, q4: theta[3],
		alpha4:  pi / 2, a4:      0, d5:     0, q5: theta[4],
		alpha5: -pi / 2, a5:      0, d6:     0, q6: theta[5],
		alpha6:     0.0, a6:      0, d7: 0.303, q7: 0}

	T0_1 = Matrix([[ cos(q1), 		  -sin(q1), 	       0,      	       a0],
				[ sin(q1 ) *cos(alpha0), cos(q1 )* cos(alpha0), -sin( alpha0), -sin(alpha0 ) *d1],
				[ sin(q1 ) *sin(alpha0) ,cos(q1) * sin(alpha0), cos(alpha0), cos(alpha0) * d1],
				[		0,			 0,	       0, 		1]])

	T0_1 = T0_1.subs(s)

	T1_2 = Matrix( [[             cos(q2),            -sin(q2),            0,              a1],
				[ sin( q2 ) *cos(alpha1), cos( q2 )* cos(alpha1), -sin(alpha1), -sin(alpha1) * d2],
				[sin(q2) * sin(alpha1), cos(q2) * sin(alpha1), cos(alpha1), cos(alpha1) * d2],
				[0, 0, 0, 1]])

	T1_2 = T1_2.subs(s)

	T2_3 = Matrix([[cos(q3), -sin(q3), 0, a2],
				[sin(q3) * cos(alpha2), cos(q3) * cos(alpha2), -sin(alpha2), -sin(alpha2) * d3],
				[sin(q3) * sin(alpha2), cos(q3) * sin(alpha2), cos(alpha2), cos(alpha2) * d3],
				[0, 0, 0, 1]])

	T2_3 = T2_3.subs(s)

	T3_4 = Matrix([[cos(q4), -sin(q4), 0, a3],
				[sin(q4) * cos(alpha3), cos(q4) * cos(alpha3), -sin(alpha3), - sin(alpha3) * d4],
				[sin(q4) * sin(alpha3), cos(q4) * sin(alpha3), cos(alpha3), cos(alpha3) * d4],
				[0, 0, 0, 1]])

	T3_4 = T3_4. subs(s)

	T4_5 = Matrix([[cos(q5), -sin(q5), 0, a4],
				[sin(q5) * cos(alpha4), cos(q5) * cos(alpha4), -sin(alpha4), -sin(alpha4) * d5],
				[sin(q5) * sin(alpha4), cos(q5) * sin(alpha4), cos(alpha4), cos(alpha4) * d5],
				[0, 0, 0, 1]])

	T4_5 = T4_5. subs(s)


	T5_6 = Matrix([[cos(q6), -sin(q6), 0, a5],
				[sin(q6) * cos(alpha5), cos(q6) * cos(alpha5), - sin(alpha5), -sin(alpha5) * d6],
				[sin(q6) * sin(alpha5), cos(q6) * sin(alpha5), cos(alpha5), cos(alpha5) * d6],
				[0, 0, 0, 1]])

	T5_6 = T5_6.subs(s)

	T6_G = Matrix([[cos(q7), -sin(q7), 0, a6],
				[sin(q7) * cos(alpha6), cos(q7) * cos(alpha6), -sin(alpha6), -sin(alpha6) * d7],
				[sin(q7) * sin(alpha6), cos(q7) * sin(alpha6), cos(alpha6), cos(alpha6) * d7],
				[0, 0, 0, 1]])

	T6_G = T6_G.subs(s)

	T0_2 = simplify(T0_1 * T1_2)
	T0_3 = simplify(T0_2 * T2_3)
	T0_4 = simplify(T0_3 * T3_4)
	T0_5 = simplify(T0_4 * T4_5)
	T0_6 = simplify(T0_5 * T5_6)
	T0_G = simplify(T0_6 * T6_G)

	R_z = Matrix([[cos(np.pi), -sin(np.pi), 0, 0],
				[sin(np.pi), cos(np.pi), 0,  0],
				[0, 0, 1  , 0],
				[      0, 	 0, 0, 1]])

	R_y = Matrix([[cos(-np.pi / 2), 0, sin(-np.pi / 2), 0],
				[	0, 1,          0 ,    0],
				[-sin(-np.pi /2), 0, cos(-np.pi /2) , 0],
				[          0, 0, 	  0, 1]])

	R_corr = simplify(R_z * R_y)

	T_total = simplify(T0_G * R_corr)

	return T_total[0,3], T_total[1,3], T_total[2,3]


def handle_calculate_IK(req):
	rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
	if len(req.poses) < 1:
		print "No valid poses received"
		return -1
	else:

		### Your FK code here
		# Create symbols
		#
		#
		# Create Modified DH parameters
		#
		#
		# Define Modified DH Transformation matrix
		#
		#
		# Create individual transformation matrices
		#
		#
		# Extract rotation matrices from the transformation matrices
		#
		#
		###
		d6_val = 0.0
		d1_val = 0.75
		a1_val = 0.35
		a2_val = 1.25
		a3_val = -0.054
		d4_val = 1.50
		l = 0.303

		q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
		d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
		a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
		alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

		s = {alpha0:    0.0, a0:     0, d1: 0.75, q1:  q1,
			alpha1:-pi / 2, a1:  0.35, d2:    0, q2:  q2-np.pi / 2,
			alpha2:    0.0, a2:  1.25, d3:    0, q3:  q3,
			alpha3:-pi / 2, a3:-0.054, d4: 1.50, q4:  q4,
			alpha4: pi / 2, a4:     0, d5:    0, q5:  q5,
			alpha5:-pi / 2, a5:     0, d6:    0, q6:  q6,
			alpha6:    0.0, a6:     0, d7:0.303, q7:  0}

		T0_1 = Matrix([[ 	cos(q1), 		  -sin(q1), 	       0,      	       a0],
				[ sin(q1)*cos(alpha0), cos(q1)* cos(alpha0), -sin( alpha0), -sin(alpha0)*d1],
				[ sin(q1)*sin(alpha0),cos(q1) * sin(alpha0), cos(alpha0), cos(alpha0) * d1],
				[		0,			 0,	       0, 		1]])

		T0_1 = T0_1.subs(s)

		T1_2 = Matrix( [[      cos(q2),            -sin(q2),            0,              a1],
				[ sin( q2)*cos(alpha1), cos( q2 )*cos(alpha1), -sin(alpha1), -sin( alpha1 )*d2],
				[ sin(q2)*sin(alpha1), cos(q2)* sin(alpha1), cos( alpha1 ),   cos(alpha1)*d2],
				[                   0,                   0,            0,               1]])

		T1_2 = T1_2.subs(s)

		T2_3 = Matrix( [[      cos(q3),            -sin(q3),            0,              a2],
				[ sin(q3)*cos( alpha2), cos(q3)*cos(alpha2), -sin(alpha2 ) , -sin(alpha2)*d3],
				[ sin(q3)*sin( alpha2), cos(q3)*sin(alpha2), cos(alpha2 ) ,   cos(alpha2)*d3],
				[                   0,                   0,            0 ,               1]])

		T2_3 = T2_3.subs(s)

		"""
		T3_4 = Matrix( [[             cos(q4),            -sin( q4),            0,              a3],
				[ sin(q4)*cos(alpha3), cos(q4)*cos(alpha3), -sin(alpha3), - sin(alpha3)*d4],
				[ sin(q4)*sin(alpha3), cos(q4)*sin(alpha3), cos( alpha3),   cos(alpha3)*d4],
				[                   0,                   0,            0,               1]])

		#T3_4 = T3_4. subs(s)

		T4_5 = Matrix( [[             cos(q5),            -sin(q5),            0,              a4] ,
				[ sin(q5)*cos(alpha4), cos( q5 )*cos(alpha4), -sin(alpha4),-sin(alpha4)*d5],
				[ sin(q5)*sin( alpha4), cos(q5)*sin(alpha4), cos(alpha4),   cos(alpha4)*d5],
				[                   0,                   0,            0 ,               1]])

		#T4_5 = T4_5. subs(s)


		T5_6 =Matrix( [[ cos( q6),            -sin(q6),            0,              a5],
				[ sin(q6)* cos(alpha5) , cos(q6)*cos(alpha5), - sin(alpha5), -sin(alpha5)*d6],
				[ sin(q6)*sin(alpha5), cos(q6)*sin( alpha5), cos( alpha5), cos( alpha5)*d6],
				[                   0,                   0,            0,               1]])

		#T5_6 = T5_6.subs(s)

		T6_G = Matrix( [[             cos(q7),            -sin(q7),            0,              a6],
				[ sin(q7)*cos(alpha6), cos(q7)*cos(alpha6), -sin(alpha6), -sin(alpha6)*d7],
				[ sin(q7)*sin(alpha6), cos(q7)*sin(alpha6),  cos(alpha6),  cos(alpha6)*d7],
				[                   0,                   0,0,               1]])

		#T6_G = T6_G.subs(s)
		"""

		T0_2 = T0_1 * T1_2
		T0_3 = T0_2 * T2_3
		#T0_4 = simplify(T0_3 * T3_4)
		#T0_5 = simplify(T0_4 * T4_5)
		#T0_6 = simplify(T0_5 * T5_6)
		#T0_G = simplify(T0_6 * T6_G)

		R_z = Matrix([[cos(np.pi), -sin(np.pi), 0],
					[  sin(np.pi),  cos(np.pi), 0],
					[      		0,	 		 0, 1]])

		R_y = Matrix([[cos(-np.pi/2), 0, sin(-np.pi/2)],
					[	   		   0, 1,          	 0],
					[ -sin(-np.pi/2), 0, cos(-np.pi/2)]])

		R_corr = simplify(R_z * R_y)

		#T_total = simplify(T0_G * R_corr)


		# Initialize service response
		joint_trajectory_list = []
		for x in xrange(0, len(req.poses)):
			# IK code starts here
			joint_trajectory_point = JointTrajectoryPoint()

			# Extract end-effector position and orientation from request
			# px,py,pz = end-effector position
			# roll, pitch, yaw = end-effector orientation

			px = req.poses[x].position.x
			py = req.poses[x].position.y
			pz = req.poses[x].position.z

			(roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
					[req.poses[x].orientation.x, req.poses[x].orientation.y,
					req.poses[x].orientation.z, req.poses[x].orientation.w])

			p1, p2, p3 = symbols('p1:4')

			Ry = rot_z(p1).subs({p1: yaw})
			Rp = rot_y(p2).subs({p2: pitch})
			Rr = rot_x(p3).subs({p3: roll})
			Rrpy = simplify(Ry * Rp * Rr * R_corr)

			wx = px - (d6_val + l) * Rrpy[0, 2]
			wy = py - (d6_val + l) * Rrpy[1, 2]
			wz = pz - (d6_val + l) * Rrpy[2, 2]

			#print wx, wy, wz

			theta1 = atan2(wy, wx)

			#print theta1

			A = sqrt(a3_val**2 + d4_val**2)
			B = sqrt((wz-d1_val)**2 + (sqrt(wy**2 + wx**2)-a1_val)**2)
			C = a2_val

			a = acos((-A ** 2 + B ** 2 + C ** 2) / (2 * C * B))
			b = acos((A ** 2 - B ** 2 + C ** 2) / (2 * A * C))
			c = acos((A**2 + B**2 - C**2)/(2*A*B))

			#print a, b, c

			theta2 = np.pi/2 - a - atan2((wz-d1_val), (sqrt(wy**2 + wx**2)-a1_val))

			theta3 = np.pi/2 - (b - atan2(a3_val, d4_val))
			#print theta1, theta2, theta3

			R0_3 = T0_3[:3, :3]
			R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

			R3_6 = R0_3.inv("LU") * Rrpy


			#R = (T3_4 * T4_5 * T5_6)[0:3, 0:3]
			#R = simplify(R)
			#print R
			"""
			for i in range(9):
				print R[i]
			    print i
			    if i+1 % 3 == 0:
			        print "\n\n\n"
			"""

			theta4 = atan2(R3_6[2, 2], -R3_6[0, 2])
			theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])

			fail = false
			theta5 = atan2(sqrt(R3_6[1, 0] ** 2 + R3_6[1, 1] ** 2), R3_6[1, 2])
			theta5_old = 0
			if np.abs(sin(theta5)*cos(theta6)-R3_6[1, 0]) > 1e-3 or np.abs(sin(theta5)*sin(theta6)+R3_6[1, 1]) > 1e-3:
				theta5_old = theta5
				theta5 = atan2(-sqrt(R3_6[1, 0] ** 2 + R3_6[1, 1] ** 2), R3_6[1, 2])
				if np.abs(sin(theta5) * cos(theta6) - R3_6[1, 0]) > 1e-3 or np.abs(
										sin(theta5) * sin(theta6) + R3_6[1, 1]) > 1e-3:
					rospy.logwarn("Solution not found for theta5")
					fail = true


			#print theta4, theta5, theta6
			### Your IK code here
			# Compensate for rotation discrepancy between DH parameters and Gazebo
			#
			#
			# Calculate joint angles using Geometric IK method
			#
			#
			###

			# Populate response for the IK request
			# In the next line replace theta1,theta2...,theta6 by your joint angle variables
			if fail:
				ee_calculated1 = forward_kinematics([theta1, theta2, theta3, theta4, theta5, theta6])
				ee_calculated2 = forward_kinematics([theta1, theta2, theta3, theta4, theta5_old, theta6])

				# print ee_calculated
				error1 = []
				error1.append(np.abs(ee_calculated1[0] - px))
				error1.append(np.abs(ee_calculated1[1] - py))
				error1.append(np.abs(ee_calculated1[2] - pz))
				overall_err1 = sqrt(error1[0] ** 2 + error1[1] ** 2 + error1[2] ** 2)

				error2 = []
				error2.append(np.abs(ee_calculated2[0] - px))
				error2.append(np.abs(ee_calculated2[1] - py))
				error2.append(np.abs(ee_calculated2[2] - pz))
				overall_err2 = sqrt(error2[0] ** 2 + error2[1] ** 2 + error2[2] ** 2)
				if overall_err1 < overall_err2:
					joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
					rospy.loginfo(str(sum(ee_calculated1)) + "\tError " + str(x) + "\t" + ": " + "\t" + str(
						error1) + "\tOverall error: " + str(overall_err1) + " other: " + str(overall_err2))

				else:
					joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5_old, theta6]
					rospy.loginfo(str(sum(ee_calculated2)) + "\tError " + str(x) + "\t" + ": " + "\t" + str(
						error2) + "\tOverall error: " + str(overall_err2)+ " other: " + str(overall_err1))

			else:
				joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
				ee_calculated = forward_kinematics(joint_trajectory_point.positions)
				#print ee_calculated
				error = []
				error.append(np.abs(ee_calculated[0] - px))
				error.append(np.abs(ee_calculated[1] - py))
				error.append(np.abs(ee_calculated[2] - pz))
				overall_err = sqrt(error[0] ** 2 + error[1] ** 2 +error[2] ** 2 )
				rospy.loginfo(str(sum(ee_calculated)) + "\tError " + str(x) + "\t" + ": " + "\t" + str(error) + "\tOverall error: " + str(overall_err))
			joint_trajectory_list.append(joint_trajectory_point)

		rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
		return CalculateIKResponse(joint_trajectory_list)


def IK_server():
	# initialize node and declare calculate_ik service
	rospy.init_node('IK_server')
	s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
	print "Ready to receive an IK request"
	rospy.spin()

if __name__ == "__main__":
	IK_server()
