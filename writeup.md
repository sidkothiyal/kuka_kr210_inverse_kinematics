## Project: Kinematics Pick & Place

[//]: # (Image References)

[image1]: ./imgs/fig.jpg
[image2]: ./imgs/img2.png
[image3]: ./imgs/img3.png

### Kinematic Analysis
#### 1. Run the forward_kinematics demo and evaluate the kr210.urdf.xacro file to perform kinematic analysis of Kuka KR210 robot and derive its DH parameters.

The figure for the Kuka KR210 was drawn, and appropriately annotated. The different joints and link properties were explored in the kr210.urdf.xacro to fill in the appropriate values in the DH parameter table.

Links | alpha(i-1) | a(i-1) | d(i) | theta(i)
---  | --- | ---  | --- | ---
0->1 | 	  0|     0| 0.75| q1
1->2 |-pi/2|  0.35|    0| -pi/2 + q2
2->3 |    0|  1.25|    0| q3
3->4 |-pi/2|-0.054| 1.50| q4
4->5 | pi/2|     0|    0| q5
5->6 |-pi/2|     0|    0| q6
6->EE| 	  0|     0|0.303| 0

![alt text][image1]

#### 2. Using the DH parameter table you derived earlier, create individual transformation matrices about each joint. In addition, also generate a generalized homogeneous transform between base_link and gripper_link using only end-effector(gripper) pose.

The homogeneous transformation matrix for each successive pair of joints can be calculated by filling in the appropriate DH parameter value, in the following matrix:

```
T(i-1)_(i) = Matrix([[ 	      cos(q(i)), 		        -sin(q(i)), 	       	   0,      	         a(i-1)],
			[ sin(q(i))*cos(alpha(i-1)), cos(q(i))*cos(alpha(i-1)), -sin(alpha(i-1)), -sin(alpha(i-1))*d(i)],
			[ sin(q(i))*sin(alpha(i-1)), cos(q(i))*sin(alpha(i-1)),  cos(alpha(i-1)),  cos(alpha(i-1))*d(i)],
			[						  0,			 			 0,	         	   0, 					  1]])
```

The above Transformation matrix for each successive joint can be post multiplied from the base_link to gripper_link to achieve the full homogeneous Transformation matrix for the gripper_link.
```
T0_G = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_G
```
If the orientation and position of the end effector is known, the homogeneous transform can also be found as follows:

```
R_corr = Rot_z(pi) * Rot_y(-pi/2)
R0_6 = Rot_z(yaw) * Rot_y(pitch) * Rot_x(roll) * R_corr

T0_G = Matrix([[R0_6[0,0], R0_6[0,1], R0_6[0,2], pos[x]],
			[	R0_6[1,0], R0_6[1,1], R0_6[1,2], pos[y]],
            [	R0_6[2,0], R0_6[2,1], R0_6[2,2], pos[z]],
            [			0,		   0,		  0, 	  1]])

```

Since the orientation of the reference frame of the base_link and gripper_link are different, the gripper needs to be transformed using intrinsic rotation about z-axis by 180 degrees and y-axis by -90 degrees.


#### 3. Decouple Inverse Kinematics problem into Inverse Position Kinematics and inverse Orientation Kinematics; doing so derive the equations to calculate all individual joint angles.

The problem of inverse kinematics aims to find the theta angle of every joint to achieve the given end effector position and orientation. In the Kuka KR210 arm, the first 3 joints are used to set the location of the wrist center (WC) and the last 3 are mostly required to define the orientation of the end effector. 
Breaking down this problem into the above two parts makes it easier to solve, since we know the end-effector position and orientation, and the length of each link can be found using the DH parameter table.

Since we know the length of the end effector and its orientation and final position, we subtract the end effector length along its axis and get the WC coordinates.

We project the WC coordinates and joint location on the x-y plane of the reference frame of joint 1, to find the angle theta1, by finding the angle the line segment from joint 1 to WC makes with x-axis of joint 1.

To compute theta2 and theta3, the joint 2, joint 3 and WC positions are projected on the x-y plane of the reference frame of joint 2. We can find the lengths of link 2 and link 3 from the DH parameter table, and the length of the line segment can be found by subtracting the position of joint 2 from the position of the WC. 

Now we have a triangle between joint 2, joint 3 and WC and we know the lengths of all the sides, and using cosine law, all three angles can be calculated. theta2 and theta3 can be found using the angles of the triangle, and small corrections need to be applied according to the DH parameter table, to keep the values consistent.

Now that we have theta1 to theta3, and the final transformation matrix (hence rotation matrix), we can compute rotation matrix R0_3, by replacing q1, q2, q3 with theta1, theta2, theta3 respectively. Using this we can find the matrix R3_6, as follows:

```
R3_6 = inv(R0_3) * R0_6
```

Simplified R3_6 looks something like this:

```
R3_6 = Matrix([[-sin(q4)sin(q6)+cos(q4)cos(q5)cos(q6),-sin(q4)cos(q6)+cos(q4)cos(q5)sin(q6),-cos(q4)sin(q5)],
			[						   sin(q5)cos(q6),						-sin(q5)sin(q6),		cos(q5)],
            [	-cos(q4)sin(q6)-sin(q4)cos(q5)cos(q6),-cos(q4)cos(q6)+sin(q4)cos(q5)sin(q6),-cos(q4)sin(q5)]])
```
Using the above matrix, theta4, theta5, theta6 are:
```
theta4 = q4 = atan2(R3_6[2,2],-R3_6[0,2])
theta5 = q5 = atan2(sqrt(R3_6[1,0]**2 + R3_6[1,1]**2), R3_6[1,2])
				or
theta5 = q5 = atan2(sqrt(R3_6[0,2]**2 + R3_6[2,2]**2), R3_6[1,2])
theta6 = q6 = atan2(-R3_6[1,1], R3_6[1,0])
```


### Project Implementation

#### 1. Fill in the `IK_server.py` file with properly commented python code for calculating Inverse Kinematics based on previously performed Kinematic Analysis. Your code must guide the robot to successfully complete 8/10 pick and place cycles. Briefly discuss the code you implemented and your results. 

The transformation matrices for the T0_3 to be computed were setup using sympy, so that values of theta can easily be changed later, and tested. Even the correlation matrix was also setup once, so that it is not initialized multiple times.

All the calculations were done with sympy and numpy, atan2 was used to calculate theta angles, so that signs are consistent. 

In case of multiple solutions, all solutions were tested and the one with minimum overall error was chosen.
In case of theta5 calculation, since sqrt of the values is being taken, the solution would be either -ve or +ve, so both options were considered.

![alt text][image2]
![alt text][image3]

