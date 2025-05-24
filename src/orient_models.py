import open3d as o3d
import numpy as np
import os
import sys
import math


'''
Derek Wingard, 03/2025

Basic idea is representing the meshes with a sparse point cloud of maximally distant points to preserve important symmetries while avoiding asymmetries which lead to inconsistent PCA. Simple orientation corrections done via ray intersection. To fine tune alignments, the meshes are projected onto a 2d eigen subspace where correction angles are found from the linear fit analysis of the subsequent 2d distributions. 
'''

CWD = os.getcwd()
input_folder = os.path.join(CWD, "raw_models")
output_folder = os.path.join(CWD,"oriented_models")




###################################################################################################


def orient_models(input_folder, output_folder):
    direction_choice = 'init'
    while direction_choice != 'head-on' and  direction_choice != 'up': 
        direction_choice = str(input("Orient with incisor facing viewer <head-on> or vertical <up>? (type an <input>) "))    
        if direction_choice != 'head-on' and  direction_choice != 'up':
            print('input 1 of 2 options')
            print('1). head-on')
            print('2). up')
        
    

    
    
    downsample_size = 2000
    i = 1

    input_folder,output_folder = str(input_folder),str(output_folder)
    

    CWD = os.getcwd()
    input_folder = os.path.join(CWD, input_folder)
    output_folder = os.path.join(CWD, output_folder)  
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    dir_len = len(os.listdir(input_folder))
    print('.')
    print('.')
    print('.')
    print('Aligning meshes . . .')
    print(' ')
    for file in os.listdir(input_folder):
        
        path_in = os.path.join(input_folder, f'sample{i}.stl')
        path_out = os.path.join(output_folder, f'sample{i}.stl')
        in_mesh = o3d.io.read_triangle_mesh(path_in)
        mesh,_= align_models(in_mesh,downsample_size)
        if direction_choice == 'up':
            theta = -np.pi/2
            axis = np.array([1,0,0])
            rotation_matrix = o3d.geometry.get_rotation_matrix_from_axis_angle(axis * theta)
            mesh.rotate(rotation_matrix, center = np.asarray(mesh.get_center()))
            theta = np.pi
            axis = np.array([0,1,0])
            rotation_matrix = o3d.geometry.get_rotation_matrix_from_axis_angle(axis * theta)
            mesh.rotate(rotation_matrix, center = np.asarray(mesh.get_center()))
            mesh.compute_vertex_normals()
        o3d.io.write_triangle_mesh(path_out, mesh)
        print(f'{i}/{dir_len}')
        print(' ')
        i+=1

     
        


    
    
        
        


###################################################################################################


def align_models(mesh,downsample_size):
    
    sample_size = 50000 
    
    mesh.compute_vertex_normals()
    vertices = np.asarray(mesh.vertices)
    point_cloud = o3d.geometry.PointCloud()        
    point_cloud.points = o3d.utility.Vector3dVector(vertices)
    num_points = len(np.asarray(point_cloud.points))
    indices = np.random.choice(num_points, sample_size, replace=False)
    sampled_points = np.asarray(point_cloud.points)[indices]
    
    furthest_partner_sampled_cloud = furthest_point_downsampler(sampled_points, downsample_size)
    sampled_point_cloud = o3d.geometry.PointCloud()

    
    
    sampled_point_cloud.points = o3d.utility.Vector3dVector(furthest_partner_sampled_cloud)
        
    sampled_point_cloud_centered, evecs =  get_principal_components( sampled_point_cloud)
        
    R = create_rotation_matrix(evecs)
    FULL_R = R
        
    sampled_point_cloud_centered_and_rotated = sampled_point_cloud_centered.rotate(R)
    _, evecs_rotated =  get_principal_components( sampled_point_cloud_centered_and_rotated)

    R = create_1D_rotation_matrix(evecs_rotated[:,2],[0,0,1])
    FULL_R = R @ FULL_R 
    sampled_point_cloud_centered_and_rotated.rotate(R)
        

    R = create_1D_rotation_matrix(evecs_rotated[:,0],[1,0,0])
    FULL_R = R @ FULL_R 
    sampled_point_cloud_centered_and_rotated.rotate(R)
    _, evecs_rotated =  get_principal_components( sampled_point_cloud_centered_and_rotated)
    R = create_1D_rotation_matrix(evecs_rotated[:,2],[0,0,1])
    FULL_R = R @ FULL_R 
    sampled_point_cloud_centered_and_rotated.rotate(R)
    _, evecs_rotated =  get_principal_components( sampled_point_cloud_centered_and_rotated)
    R = create_1D_rotation_matrix(evecs_rotated[:,1],[0,0,1])
    FULL_R = R @ FULL_R 
    sampled_point_cloud_centered_and_rotated.rotate(R)
    _, evecs_rotated =  get_principal_components( sampled_point_cloud_centered_and_rotated)
        
    
    
    mesh.vertices=o3d.utility.Vector3dVector(np.asarray(mesh.vertices)-np.asarray(mesh.get_center()))
    FULL_R /= np.linalg.norm(FULL_R)
    mesh.rotate(FULL_R)

    check = check_ray_intersects_mesh(mesh,[0,0,1])

    if check == False :
        R = o3d.geometry.get_rotation_matrix_from_axis_angle(np.array([0,1,0]) * np.pi)
        mesh.rotate(R)
        sampled_point_cloud_centered_and_rotated.rotate(R)

    check = check_ray_intersects_mesh(mesh,[0,-1,0])

    if check == True :
        R = o3d.geometry.get_rotation_matrix_from_axis_angle(np.array([0,0,1]) * np.pi)
        mesh.rotate(R)
        sampled_point_cloud_centered_and_rotated.rotate(R)

    mesh.compute_triangle_normals
    triangle_normals = np.asarray(mesh.triangle_normals)
    check = np.sum(triangle_normals @ np.array([0,1,0]))

    if check > 0 :
        R = o3d.geometry.get_rotation_matrix_from_axis_angle(np.array([0,0,1]) * np.pi)
        mesh.rotate(R)
        sampled_point_cloud_centered_and_rotated.rotate(R)

        

    copy_mesh = o3d.geometry.TriangleMesh() 
    copy_mesh.vertices = mesh.vertices
    copy_mesh.triangles = mesh.triangles
    copy_mesh.vertex_normals = mesh.vertex_normals
    copy_mesh.triangle_normals = mesh.triangle_normals

    _,_,slope = flatten_component(copy_mesh,2)
    theta = -np.arctan(slope) 
        
    R = o3d.geometry.get_rotation_matrix_from_axis_angle(np.array([0,0,1]) * theta)
    mesh.rotate(R)
    sampled_point_cloud_centered_and_rotated.rotate(R)

        

    copy_mesh = o3d.geometry.TriangleMesh() 
    copy_mesh.vertices = mesh.vertices
    copy_mesh.triangles = mesh.triangles
    copy_mesh.vertex_normals = mesh.vertex_normals
    copy_mesh.triangle_normals = mesh.triangle_normals

    R,slope, sampless, _ = correct_forward_tilt_angle_matrix(mesh,20)
    
    
    
    mesh.rotate(R)
    sampled_point_cloud_centered_and_rotated.rotate(R)
       
    _, evecs_rotated =  get_principal_components( sampled_point_cloud_centered_and_rotated)
        
    mesh.compute_vertex_normals()

    test_x = np.abs(np.array([1,0,0]) @ evecs_rotated[:,0])
    test_y = np.abs(np.array([0,1,0]) @ evecs_rotated[:,2])
    test_z = np.abs(np.array([0,0,1]) @ evecs_rotated[:,1])

    err = 1-(test_x+test_y+test_z)/3
    
       
    return [mesh, err]
    
def get_principal_components(point_cloud):
    
    points = np.asarray(point_cloud.points)
    center = np.mean(points, axis=0)
    
    points_centered =  points - center
    
    
    cov_matrix = np.cov(points_centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
    idx = eigenvalues.argsort()[::-1]
    eigenvectors = eigenvectors[:, idx]
    eigenvectors = orient_pca_consistently(eigenvectors)
  
    centered_point_cloud = o3d.geometry.PointCloud()
    centered_point_cloud.points = o3d.utility.Vector3dVector(points_centered)
    
    return [centered_point_cloud, eigenvectors]
    
def orient_pca_consistently(eigenvectors):
    unit = [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])]
    for i in range(3):
        proj = np.dot(unit[i], eigenvectors[:, i])
        if proj < 0:
            eigenvectors[:, i] = -eigenvectors[:, i]
    return eigenvectors

def create_rotation_matrix(eigenvectors):

    rot_mat = np.c_[eigenvectors[:, 0],eigenvectors[:, 1],eigenvectors[:, 2]]
    rot_mat /= np.linalg.norm(rot_mat)
    return rot_mat

def create_1D_rotation_matrix(eigenvector,axis):

    axis = np.asarray(axis)
    eigenvector /= np.linalg.norm(eigenvector)
    rot_axis = np.cross(eigenvector,axis)
    theta = np.arccos(axis @ eigenvector )
    
    rot_mat = o3d.geometry.get_rotation_matrix_from_axis_angle(rot_axis * theta)
    
    
    
    rot_mat /= np.linalg.norm(rot_mat)
    return rot_mat

def check_ray_intersects_mesh(mesh, vector):

    origin = np.array([0,0,0])
    direction = np.array(vector)
    
    mesh_for_ray = o3d.t.geometry.TriangleMesh.from_legacy(mesh)
    scene = o3d.t.geometry.RaycastingScene()
    scene.add_triangles(mesh_for_ray)

    rays = o3d.core.Tensor([[*origin, *direction]], dtype=o3d.core.Dtype.Float32)
    intersection_counts = scene.count_intersections(rays).numpy()
    
    return intersection_counts[0] > 0


def flatten_component(mesh,comp):
    
    R = np.eye(3)

    R[comp,:] *= 0

    mesh.rotate(R)

    av,slope = create_alignment_vector(mesh,comp)
    
    return [mesh,av,slope]



def create_alignment_vector(mesh,comp):
    idx1, idx2 = (comp+1)%3, (comp+2)%3

    e1,e2 = np.asarray(mesh.vertices)[:,idx1],np.asarray(mesh.vertices)[:,idx2]
    sample_size = int(len(e1)*0.001)
    indices = np.random.choice(len(e1), size=sample_size, replace=False)
    e1_data, e2_data = e1[indices], e2[indices]

    coefficients = np.polyfit(e1_data, e2_data, 1)
    
    a = coefficients[0]
  
    alignment_vec = np.array([0,0,0])
    alignment_vec[idx1] = a
    alignment_vec[idx2] = 1
    

    return [alignment_vec, a]

def correct_forward_tilt_angle_matrix(mesh,n):
    
    vertices = np.asarray(mesh.vertices)
    z = vertices[:,2]
    sorted_indices = np.argsort(z)
    total = len(vertices)
    lower_thresh = 1000
    upper_thresh = total - lower_thresh
    lower_indices = sorted_indices[:lower_thresh]
    upper_indices = sorted_indices[upper_thresh:]

    n2 = int(n/2)
    
    lower_sampled_indices = np.random.choice(lower_indices, n2 , replace=False)
    upper_sampled_indices = np.random.choice(upper_indices, n2 , replace=False)

    lower_samples = vertices[lower_sampled_indices]
    upper_samples = vertices[upper_sampled_indices]

    lower_samples[:,0] *= 0
    upper_samples[:,0] *= 0

    lower_average = np.mean(lower_samples, axis=0)
    upper_average = np.mean(upper_samples, axis=0)

    samples = np.concatenate([lower_samples, upper_samples])

    

    vector = upper_average - lower_average 

    vector /= np.linalg.norm(vector)

    axis = np.cross(vector,np.array([0,0,1]) )
    theta = np.arccos(vector @ np.array([0,0,1]))

    slope = vector[1]/vector[2]
    theta = 0.5*np.arccos(slope)

    R = o3d.geometry.get_rotation_matrix_from_axis_angle(axis * theta)

    return [R,slope,samples,vector]

    

    
    
def furthest_point_downsampler(point_cloud, num_points):

    points = np.asarray(point_cloud)
    num_samples = min(num_points, points.shape[0])
    sampled_indices = np.zeros(num_samples, dtype=np.int32)
    sampled_indices[0] = np.random.randint(0, points.shape[0])
    used_points = np.zeros(points.shape[0], dtype=bool)
    start_idx = np.random.randint(0, points.shape[0])
    sampled_indices[0] = start_idx
    used_points[start_idx] = True
    for i in range(1, num_samples):
        # Get the last sampled point
        last_point = points[sampled_indices[i-1]]
        
        # Calculate squared distances from all points to the last sampled point
        distances = np.sum((points - last_point)**2, axis=1)
        distances[used_points] = -1 
        
        # Find the furthest point
        furthest_point_idx = np.argmax(distances)
        
        if distances[furthest_point_idx] < 0:
            print(f"Warning: Only able to sample {i} unique points")
            num_samples = i
            break
        
        # Store the index of the furthest point
        sampled_indices[i] = furthest_point_idx
        
        used_points[furthest_point_idx] = True
        
    sampled_indices = sampled_indices[:num_samples]
    furthest_point_downsampled_cloud = point_cloud[sampled_indices]
    
    return furthest_point_downsampled_cloud


orient_models(input_folder, output_folder)

rand = np.random.randint(1,11)
oriented_file = os.path.join(output_folder, f'sample{rand}.stl')
mesh = o3d.io.read_triangle_mesh(oriented_file)
mesh.compute_vertex_normals()
o3d.visualization.draw_geometries([mesh])




