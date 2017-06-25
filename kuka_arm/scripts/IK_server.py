#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2
from sympy.matrices import Matrix
import pickle, os


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []

        def DH_transform(q1, d1, alpha0, a0): 
            '''The modified DH convention transform matrix
            alpha0: twist angle, a0: link length, 
            d1: link offset, q1: joint angle'''
            T0_1 = Matrix([[cos(q1), -sin(q1),      0,    a0],
                           [sin(q1) * cos(alpha0), cos(q1) * cos(alpha0), -sin(alpha0), -sin(alpha0) * d1],
                           [sin(q1) * sin(alpha0), cos(q1) * sin(alpha0),  cos(alpha0),  cos(alpha0) * d1],
                           [0,              0,      0,     1]])
            return T0_1

        # def rot_x(q):
        #     M1 = Matrix([[0],[0],[0]])
        #     M2 = Matrix([[0, 0, 0, 1]])
        #     Rx = Matrix([[  1,      0,       0],
        #                  [  0, cos(q), -sin(q)],
        #                  [  0, sin(q),  cos(q)]])
        #     rot_x = Rx.row_join(M1)
        #     rot_x = rot_x.col_join(M2)
        #     return rot_x
            
        # def rot_y(q):   
        #     M1 = Matrix([[0],[0],[0]])
        #     M2 = Matrix([[0, 0, 0, 1]])
        #     Ry = Matrix([[ cos(q),     0, sin(q)],
        #                  [      0,     1,      0],
        #                  [-sin(q),     0, cos(q)]])
        #     rot_y = Ry.row_join(M1)
        #     rot_y = rot_y.col_join(M2)
        #     return rot_y

        # def rot_z(q):   
        #     M1 = Matrix([[0],[0],[0]])
        #     M2 = Matrix([[0, 0, 0, 1]])
        #     Rz = Matrix([[ cos(q),-sin(q), 0],
        #                  [ sin(q), cos(q), 0],
        #                  [      0,      0, 1]])
        #     rot_z = Rz.row_join(M1)
        #     rot_z = rot_z.col_join(M2)
        #     return rot_z
        def rot_x(q):
            R_x = Matrix([[      1,      0,      0],
                          [      0, cos(q), -sin(q)],
                          [      0, sin(q),  cos(q)]])
            
            return R_x
    
        def rot_y(q):
            R_y = Matrix([[ cos(q),     0, sin(q)],
                          [      0,     1,      0],
                          [-sin(q),     0, cos(q)]])
            
            return R_y

        def rot_z(q):
            R_z = Matrix([[ cos(q),-sin(q), 0],
                          [ sin(q), cos(q), 0],
                          [      0,      0, 1]])
            
            return R_z

        def Euler_angles_from_rotation_matrix(R):
            '''Input R is 3x3 matrix, output are alpha, beta and gamma
            alpha: rotation angle along z axis,
            beta: rotation angle along y axis
            gamma: rotation angle along x axis'''

            r11, r21, r31 = R[0,0], R[1,0], R[2,0]
            r32, r33 = R[2,1], R[2,2]

            # Euler angles from rotation matrix
            beta = atan2(-r31, sqrt(r11**2 + r21**2))
            gamma = atan2(r32, r33)
            alpha = atan2(r21, r11)
            return alpha, beta, gamma


        theta1,theta2,theta3,theta4,theta5,theta6 = 0,0,0,0,0,0
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
        r, p, y = symbols("r p y") # end-effector orientation
        R0_g_symp = simplify(rot_x(r) * rot_y(p) * rot_z(y))

        # The Modified DH params
        s = {alpha0:       0, a0:      0, d1:  0.75,
             alpha1: -pi / 2, a1:   0.35, d2:     0, q2: q2 - pi / 2,
             alpha2:       0, a2:   1.25, d3:     0,
             alpha3: -pi / 2, a3: -0.054, d4:   1.5,
             alpha4:  pi / 2, a4:      0, d5:     0,
             alpha5: -pi / 2, a5:      0, d6:     0,
             alpha6:       0, a6:      0, d7: 0.303, q7: 0}

        # Define Modified DH Transformation matrix
        if not os.path.exists("R0_3_inv.p"):
            print("-----------------------------------------")
            print("No DH matrices, create and save it.")
            # base_link to link1
            T0_1 = DH_transform(q1, d1, alpha0, a0)
            T0_1 = T0_1.subs(s)
            # linke1 to link 2
            T1_2 = DH_transform(q2, d2, alpha1, a1)
            T1_2 = T1_2.subs(s)
            # link2 to link3
            T2_3 = DH_transform(q3, d3, alpha2, a2)
            T2_3 = T2_3.subs(s)
            # link3 to link4
            # T3_4 = DH_transform(q4, d4, alpha3, a3)
            # T3_4 = T3_4.subs(s)
            # # link4 to link5
            # T4_5 = DH_transform(q5, d5, alpha4, a4)
            # T4_5 = T4_5.subs(s)
            # # link5 to link6
            # T5_6 = DH_transform(q6, d6, alpha5, a5)
            # T5_6 = T5_6.subs(s)
            # # link6 to end-effector
            # T6_7 = DH_transform(q7, d7, alpha6, a6)
            # T6_7 = T6_7.subs(s)

            # Create individual transformation matrices
            T0_2 = simplify(T0_1 * T1_2)
            # R0_3 inv matrix would be used in R3_6
            T0_3 = simplify(T0_2 * T2_3)
            R0_3 = T0_3[0:3, 0:3]
            R0_3_inv = simplify(R0_3 ** -1)
            # T0_4 = simplify(T0_3 * T3_4)
            # T0_5 = simplify(T0_4 * T4_5)
            # T0_6 = simplify(T0_5 * T5_6)
            # T0_G = simplify(T0_6 * T6_7)

            # R_corr = rot_z(pi) * rot_y(-pi/2)
            #T_total = simplify(T0_G * R_corr)
            pickle.dump(T0_2, open("T0_2.p", "wb"))
            pickle.dump(R0_3_inv, open("R0_3_inv.p", "wb"))
            print("Transformation matrices have been saved!!")
            print("-----------------------------------------")
        else:
            T0_2 = pickle.load(open("T0_2.p", "rb"))
            R0_3_inv = pickle.load(open("R0_3_inv.p", "rb"))
            print("-----------------------------------------")
            print("Transformation matrices have been loaded!")
            print("-----------------------------------------")           
        
        # Compensation matix for wrist center calculation
        R_corr = rot_z(pi) * rot_y(-pi/2)
        print("Calculating Inverse kinematrics")

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
     
            # Calculate joint angles using Geometric IK method
            pg_0 = Matrix([[px],[py],[pz]])
            R0_g_value = R0_g_symp.evalf(subs={r: roll, p: pitch, y: yaw})
            R0_g = R0_g_value[0:3,0:3]

            # Calculate wrist center and theta1
            pwc_0 = pg_0 - (0.303) * R0_g * R_corr * Matrix([[0],[0],[1]])
            # pwc_0 = pg_0 - (0.303) * R0_g * Matrix([[0],[0],[1]])
            theta1 = atan2(pwc_0[1], pwc_0[0])

            # Calculate theta3
            p2_sym = T0_2 * Matrix([0,0,0,1])
            p2_0 = p2_sym.evalf(subs={q1: theta1})
            pwc_2 = pwc_0 - p2_0[0:3,:]

            l23 = a2
            l35 = sqrt(a3**2 + d4**2)
            p25 = sqrt(np.sum(np.square(pwc_2)))
            theta3_1 = atan2(a3,d4)
            c235 = (np.sum(np.square(pwc_2)) - l23**2 - l35**2) / (2*l23*l35)
            theta3_phi = atan2(sqrt(1-c235**2), c235)
            # theta3_phi = atan2(-sqrt(1-c235**2), c235)
            theta3 = (theta3_phi + theta3_1 - pi/2).subs(s)

            # Calculate theta2
            theta2_out = atan2(pwc_2[2], sqrt(pwc_2[0]**2 + pwc_2[1]**2))
            c523 = (-l35**2 + l23**2 + p25**2) / (2*l23 * p25)
            theta2_phi = atan2(sqrt(1 - c523**2), c523)
            theta2 = (pi/2 - (theta2_phi + theta2_out)).subs(s)

            # Calculate theta4, 5, 6
            R0_3_inv_value = R0_3_inv.evalf(subs={q1:theta1, q2:theta2, q3:theta3})
            R3_6 =  R0_3_inv_value * R0_g

            alpha, beta, gamma = Euler_angles_from_rotation_matrix(R3_6)
            theta4 = gamma
            theta5 = beta
            theta6 = alpha
        
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)
        print("Inverse kinematics calculation has been done!")
        print("-----------------------------------------")
        rospy.loginfo("\nlength of Joint Trajectory List: %s" % len(joint_trajectory_list))
        # rospy.loginfo("Inverse kinematics calculation has been done!")
        # rospy.loginfo("Px, Py, Pz %s %s %s" , px , py , pz)
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
