# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:59:38 2024

@author: Admin
"""
import numpy as np
import os

L=64
#L=32
N=L**2
restlen=np.ones(shape=((L**2),6))
p=1
bond_len = 1
dipole_number=5
#r=6
r=25
#center=120
#center=496
center=2080
final_restlen=0.9
total_step=10000
#total_step=1    
#the connection matrix for each node
def connect_matrix(L):
    connect=-1*np.ones((L**2,6), dtype=int)
    n=0
    for n in range(L**2):
        connect[n][0]=n-1
        connect[n][1]=n-L-1
        connect[n][2]=n-L
        connect[n][3]=n+1
        connect[n][4]=n+L
        connect[n][5]=n+L-1
    #left edges even rows
        if n%L==0:
            connect[n][0]=-1
    #even rows        
        if (int(n/L)%2)==0:
            connect[n][1]=n-L
            connect[n][2]=n-L+1
            connect[n][4]=n+L+1
            connect[n][5]=n+L
            #bottom nodes      
            if n<L:
                connect[n][1]=-1
                connect[n][2]=-1
            #top nodes
        if (n+L) >= (L**2):
            connect[n][4]=-1
            connect[n][5]=-1
    #left edges, odd rows        
        if (((n/L)-1))%2==0:
            connect[n][1]=-1
            connect[n][5]=-1
    #right edges odd rows
        if (n+1)%L==0:
            connect[n][3]=-1
    #right edges, even rows
        if (((n+1)/L)-1)%2==0:
            connect[n][2]=-1
            connect[n][4]=-1
    #connect_full=connect.copy()
    return  connect
connect=connect_matrix(L)  
#************************************************************************************************************
def create_dipoles(N): 
    #rand_dipole = random.randint(min_val1, max_val1)
    dip_center =[center]
    dip_node=[center]
    for n in range (0,N):
        for m in range (0,6):
            for c in range(0,len(dip_center)):
                if n==dip_center[c]:
                        dip_node.append(connect[n][m])
    return dip_center, dip_node
dip_center, dip_node=create_dipoles(N)   
#************************************************************************************************************
 # initialize the network at some positions xcord and ycord
def initial_network(L):
    xcord=np.zeros(shape=(L**2))
    for i in range(N):
        #for r in range(L):
            if (i//L)%2==0:
                xcord[i] =i%L+bond_len*1.5
            else:
                xcord[i]= i%L+bond_len
    #j=[]
    initial_value = 1
    ycord=np.zeros(shape=(L**2))
    for j in range(N):
        if j<L:
            ycord[j]=initial_value
        else:
            ycord[j] = initial_value +(j//L)*0.5*np.sqrt(3)
    return xcord, ycord

xcord,ycord=initial_network(L)


#***********************************************************************************
def circle_boundary(xcord,ycord):
    circle_boundary_node=[]
    inner_circle=[]
    connect_circle=connect.copy()
    radius=[]
    #outer_circle=[]
    x_radius=xcord-xcord[center]
    y_radius=ycord-ycord[center]
    radius=(np.sqrt(x_radius**2 + y_radius**2))
    for i in range(N):
            if radius[i]<=r-0.5:
              inner_circle.append(i)  
            #if (r-1)<=radius<=r:
            if (r-0.5)<=radius[i]<=r+0.5:
              circle_boundary_node.append(i)
            for j in range(6):
                    if radius[connect_circle[i][j]]>r+0.5 or radius[i]>r+0.5:
                        #outer_circle.append(i)
                        connect_circle[i][j]=-1
    return radius,inner_circle, circle_boundary_node,connect_circle
radius,inner_circle, circle_boundary_node,connect_circle=circle_boundary(xcord,ycord)
def  restlen_change(restlen):
    # final_restlen=0.8
    for n in range (0,L**2):
        for m in range (0,6):
            for c in range(0,len(dip_center)):
                if connect_circle[n][m]==dip_center[c]:
                    restlen[n][m]=final_restlen
    #Drestlen[dip_center,:]=0.01 
    restlen[dip_center,:]=final_restlen      
    #if restlen[dip_center[0]][0]>final_restlen:
     #   restlen-=Drestlen
    return restlen
restlen=restlen_change(restlen)
base_directory = f'E:/FIRE_algorithm/Test_run2/max steps={total_step}_test_circle{L}/{dipole_number}dipoles/'
input_directory = f'E:/FIRE_algorithm/Test_run2/max steps={total_step}_test_circle{L}/{dipole_number}dipoles/'

fx_file = input_directory + 'Fx.npy'
fy_file = input_directory + 'Fy.npy'

xpos_file = input_directory + 'xnew_.npy'
ypos_file = input_directory + 'ynew_.npy'

fx_data = np.loadtxt(fx_file)
fy_data = np.loadtxt(fy_file)
xpos_data = np.loadtxt(xpos_file)
ypos_data = np.loadtxt(ypos_file)

fx_boundary = fx_data[circle_boundary_node]
fy_boundary = fy_data[circle_boundary_node]
xpos_boundary = xpos_data[circle_boundary_node]
ypos_boundary = ypos_data[circle_boundary_node]

center_x = xpos_data[center]
center_y = ypos_data[center]

radial_vec = []
for i in range(0,np.shape(xpos_boundary)[0]):
    radial_x = center_x - xpos_boundary[i]
    radial_y = center_y - ypos_boundary[i]
        
    radial_vec.append([radial_x,radial_y])

radial_vec = np.array(radial_vec)

bndry_force_vec = np.column_stack((fx_boundary, fy_boundary))
xx_moment_sum_array = bndry_force_vec[:,0] * radial_vec[:,0]    # fx dot rx
yy_moment_sum_array = bndry_force_vec[:,1] * radial_vec[:,1]    # fy dot ry

xx_moment_sum = np.sum(xx_moment_sum_array)
yy_moment_sum = np.sum(yy_moment_sum_array)
dfar_fire = xx_moment_sum + yy_moment_sum     # this is the dipole moment at the boundary

radial_force=(xx_moment_sum_array+yy_moment_sum_array)/(np.sqrt(radial_vec[:,0]**2+ radial_vec[:,1]**2))
sum_rF=sum(radial_force)
all_bnd_force=np.column_stack((circle_boundary_node,fx_boundary,fy_boundary,radial_force))

file_bnd_force=f"boundary_forces{total_step}.npy"
total_bndF=f"boundary_totalF{total_step}.npy"
file_path_bnd_force = os.path.join(base_directory, file_bnd_force)
file_path_total_bndF= os.path.join(base_directory, total_bndF)
heading = '# nodes       Fx              Fy              Radial F'
fmt = '%12d', '%15.7e', '%15.7e', '%15.7e'
np.savetxt(file_path_bnd_force,all_bnd_force, header = heading, fmt = fmt)
np.savetxt(file_path_total_bndF,[sum_rF])
heading_po='# nodes       x              y              Radial F'
heading = '# nodes       Fx              Fy              Radial F'
fmt = '%12d', '%15.7e', '%15.7e', '%15.7e'
#np.savetxt(force_outfname, np.column_stack((boundary_nodes, force_boundary_x[-1], force_boundary_y[-1], perp_boundary_force[-1])), header = heading, fmt = fmt)
#**********************************************************************************

# d_loc calculations:: it is the local dipole moment

import math
mu = 1
mu_c = 1
def angle(x1,y1,x2,y2):
    dx = x2-x1
    dy = y2-y1
    return(math.atan2(dy,dx))

# validating that the above short-cut is correct by finding unbalanced forces on the outer dipole nodes
y_force_unbalanced = []
x_force_unbalanced = []
strain_loc = []
test_force_list = []
dist_dip_node_list = []
bonds_list = [0,1,2,3,4,5]

j = 0
for i in range(0,len(dip_node)):
    if i%7 == 0: 
        last_center = dip_node[i]
        index = int(i/7)
        y_force_unbalanced.append([])
        x_force_unbalanced.append([])
        dip_nodes_outer = dip_node[i+1:i+7]
        test_force_list.append([])
        
        
        for bond in bonds_list:
            conn_val = connect[dip_node[i]][bond]     # node that is the fourth bond connnected to
            x1 = xpos_data[dip_node[i]]
            y1 = ypos_data[dip_node[i]]
            x2 = xpos_data[conn_val]
            y2 = ypos_data[conn_val]
            dx = x1-x2
            dy = y1-y2
            dist = np.sqrt(dx**2+dy**2) 
            dist_dip_node_list.append(dist)
            print('dist is: ', dist)

                
    if i%7 != 0:

        # find the orientation of the spring connected to dipole node
        # x_last_center = xpos[-1,last_center]
        # y_last_center = ypos[-1,last_center]
        x_last_center = center_x
        y_last_center = center_y
        # x_node = xpos[-1,dip_nodes[i]]
        # y_node = ypos[-1,dip_nodes[i]]
        x_node = xpos_data[dip_node[i]]
        y_node = ypos_data[dip_node[i]]
        
        spring_vec = [x_node-x_last_center, y_node-y_last_center]
        
        angle_spring = angle(x_node, y_node, x_last_center, y_last_center)
        
        force_node_y = 0
        force_node_x = 0
        test_force_node = 0
        
        for bond in bonds_list:
            # strain_val = strain[-1][dip_node[i]][bond]    # fourth bond strain
            # conn_val = conn_node[dip_node[i]][bond]     # node that is the fourth bond connnected to
            conn_val = connect[dip_node[i]][bond]     # node that is the fourth bond connnected to
            
            # if np.isnan(strain_val) == False:# and conn_val != last_center:
            if conn_val != -1: # and conn_val != last_center:
                
                # print(dip_nodes[i], conn_val, strain_active_force)
                
                x1 = xpos_data[dip_node[i]]
                y1 = ypos_data[dip_node[i]]
                x2 = xpos_data[conn_val]
                y2 = ypos_data[conn_val]
                dx = x1-x2
                dy = y1-y2
                dist = np.sqrt(dx**2+dy**2) 
                strain_val = dist - restlen[dip_node[i],bond]
                #print('here $$',dip_node[i], bond, dist,restlen[dip_node[i],bond])
                angle_val = angle(x1, y1, x2, y2)

            if conn_val == last_center:
                # strain_active_force = 1. - (rlen + strain_val)
                strain_active_force = 1. - (final_restlen + strain_val)
                # print(dip_nodes[i], bond, strain_active_force, strain_val)
                print('here[]',dip_node[i], bond, strain_active_force,final_restlen, strain_val)
            else:
                strain_active_force = strain_val
                
            #print(dip_node[i], bond, strain_active_force)
            
            if strain_val >= 0:
#                    force_bond = -1*strain_val*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
#                    force_bond = strain_val*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
                force_bond = strain_active_force*mu                 # multiplying by -1 so that tensile bonds connected to top and bottom edges are applying negative force (contraction) on the boundary
            elif strain_val < 0:
#                    force_bond = -1*strain_val*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary
#                    force_bond = strain_val*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary
                force_bond = strain_active_force*mu_c               # multiplying by -1 so that compressive bonds connected to top and bottom edges are applying positive force (expansion) on the boundary

            force_bond_y = force_bond*np.sin(angle_val)
            force_bond_x = force_bond*np.cos(angle_val)
            
            # if dip_nodes[i] == 2079:
            #     print(bond, force_bond_x, force_bond_y, np.sqrt(np.square(force_bond_x) + np.square(force_bond_y)), force_bond)
            
            force_node_y = force_node_y + force_bond_y  # adding bond forces to force on nodes
            force_node_x = force_node_x + force_bond_x  # adding bond forces to force on nodes
            
            force_unbalanced_temp = force_bond*np.cos(angle_spring)
            # finding the projection of spring force along the radial outward vector                
            test_force = np.dot([force_bond_x, force_bond_y], spring_vec) / np.linalg.norm(spring_vec)
            # test_force = np.dot([force_bond_x, force_bond_y], spring_vec)
            # test_force = 
            
            if i%7 == 1 and bond == 3:
                test_force = -1*test_force
            if i%7 == 2 and bond == 4:
                test_force = -1*test_force
            if i%7 == 3 and bond == 5:
                test_force = -1*test_force
            if i%7 == 4 and bond == 0:
                test_force = -1*test_force
            if i%7 == 5 and bond == 1:
                test_force = -1*test_force
            if i%7 == 6 and bond == 2:
                test_force = -1*test_force
            
            if dip_node[i] == 2079:
                print([force_bond_x, force_bond_y])#, spring_vec)
                # print(bond, force_bond, angle_val, force_unbalanced_temp, test_force, np.dot([force_bond_x, force_bond_y], spring_vec), np.linalg.norm(spring_vec))
                # print(test_force)
            
            test_force_node = test_force_node + test_force
                    
        test_force_list[index].append(test_force_node)
                                
        # y_force_unbalanced[index].append(force_node_y)
        # x_force_unbalanced[index].append(force_node_x)
        
#        y_force_unbalanced[i].append(force_node_y)
#        x_force_unbalanced[i].append(force_node_x)

        j += 1

test_force_list = np.array(test_force_list)

# y_force_unbalanced = np.array(y_force_unbalanced)
# x_force_unbalanced = np.array(x_force_unbalanced)

# force_mag_unbalanced_list = np.sqrt(np.square(x_force_unbalanced) + np.square(y_force_unbalanced))
# force_mag_unbalanced = np.sum(force_mag_unbalanced_list)

# finding dipole moment using unblanced forces on the outer nodes
# d_loc_unbalanced_list = dist_dip_node_list*force_mag_unbalanced_list
# d_loc_unbalanced = np.sum(d_loc_unbalanced_list)

force_list=np.reshape(test_force_list,-1)
d_loc_test_array = dist_dip_node_list*force_list
d_loc_test = np.sum(d_loc_test_array)

print('D_loc is : ', d_loc_test)
print('D_far is : ', dfar_fire)
file_dfar= "dfar.npy"
file_d_loc_test = "d_loc_test.npy"

file_path_dfar = os.path.join(base_directory, file_dfar)
file_path_d_loc_test = os.path.join(base_directory, file_d_loc_test)

# %%
np.savetxt(file_path_dfar, [dfar_fire])

np.savetxt(file_path_d_loc_test, [d_loc_test])
