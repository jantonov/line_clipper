
#Line Clipper Version 1.0 October 15, 2016
#Author Jordan A. Sahattchieve

# Basics: throughout this code a point in 2D/3D is a pair/tripple of float numbers.
# A line segment is an ordered pair of 2D/3D points.
# The function clipping() will take as arguments a line segment l, an ordered list of points F and a float number h.  The list 
# F should represent the vertices of a (not-necessarily convex) polygon in 3D.  clipping() will excise the parts of the image of l under the perspective projection through C=(0,0,0) on the view plane
# z=h which are hidden from view w.r.t. C due to obstruction of visibility caused by F, and will return the list of segments which correspond to the unoccluded parts of l.
# The function will not clip l if F is not contained in a unique affine plane in 3D (i.e. if F is colinear or a point, or if F spans a 3-dimensional polytope).
#CAUTION: Make sure that l lies above the image plane and that l does not intersect F!  Intersections with the affine plane supporting F is OK and is handled by the current version of the function.
#CAUTION: The case when the perspective projection of l intersects the perspective projection of F nontransversally is not handled. This is left for version 2.0.
#NB: This code has not been optimized for either speed or memory and was written as a coding exercise to be embedded into a larger project of Geomagical Labs.  The author is not concerned with the code's layout (as long as it is readable), but will greatly appreciate comments on its functionality, and the algorithm's correctness!


import numpy as np
import math

def clipping(l,F,h):

    #intersection() will return the point of intersection if there is an intersection between the line segments seg_0 and seg_1 in 2D and None otherwise. seg_0 is assumed to be disjoint or transversal to seg_1.
    def intersection(seg_0,seg_1):
	s0_0=np.asarray(seg_0[0])
        s0_1=np.asarray(seg_0[1])
        s1_0=np.asarray(seg_1[0])
        s1_1=np.asarray(seg_1[1])
        t_num=np.transpose(np.array([s0_0-s1_0,s0_0-s0_1]))
	s_num=np.transpose(np.array([s1_1-s1_0,s0_0-s1_0]))
        den=np.transpose(np.array([s1_1-s1_0,s0_0-s0_1]))
        den_det=np.linalg.det(den)
        if (den_det==0): return None
        t=np.linalg.det(t_num)/den_det
        s=np.linalg.det(s_num)/den_det
        if (((t>=0) and (t<=1)) and ((s>=0) and (s<=1))):
           return (seg_1[0][0]+t*(seg_1[1][0]-seg_1[0][0]),seg_1[0][1]+t*(seg_1[1][1]-seg_1[0][1]))
	else: None 
        
    #Returns the dimension of the smallest affine subspace which contains all the points in the list F 
    def dim_F(F):
	F_vec=[]
	for i in range(1,len(F)):
	    F_vec.append(np.asarray(F[i])-np.asarray(F[0]))
	F_vec_mat=np.matrix(F_vec)
	return np.linalg.matrix_rank(F_vec_mat)
    #Break and exit if F is not a 2-dimensional polygon in 3D.
    if (dim_F(F)!=2):
	print 'Attention: The list to vertices in 3D does not span a proper 2-dimensional polygon'
	return l       
    #Implements perspective projection with center (0,0,0) and image plane z=h in 3D; the optical axis is the z-axis.  Arguments are a list L of points in 3D to be projected and the distance h of the image plane from xy-plane. The function returns a list of projected 2D points in the image plane.  
    def proj(L,h):
	pr_L=[]	
	for p in L:
	    if (p[2]<=h):
	       print "Some of the points used in proj() are behind the image plane!  proj() will return an empty list and stop."
	       return []
	for p in L:
	    pr_L.append((p[0]*h/p[2],p[1]*h/p[2]))
	return pr_L
    #Implements an inversion of the perspective projection on lines: given a line segment seg in 3D and a point P in its projected image in 2D, it returns the xyz coordinates of the preimage of the point P in seg.
    def preim(seg,P):
	a=np.asarray(seg[0])
	b=np.asarray(seg[1])
	p=np.asarray((P[0],P[1],h))
	abp_cross=np.cross(b-a,p)
	if (abp_cross[0]==0 and abp_cross[1]==0 and abp_cross[2]==0):
	   print "Error: Attempting to call preim() with a segment and a point not in its perspective projection image."
	   return None
	if (abp_cross[2]!=0):
	   den=(b-a)[1]*p[0]-(b-a)[0]*p[1]
	   num=a[0]*p[1]-a[1]*p[0]
	   t=num/den
	   if ((t>=0) and (t<=1)):
	      pre_img=a+t*(b-a)
	      return (pre_img[0],pre_img[1],pre_img[2])
	if (abp_cross[1]!=0):
	   den=(b-a)[2]*p[0]-(b-a)[0]*p[2]
	   num=a[0]*p[2]-a[2]*p[0]
	   t=num/den
	   if ((t>=0) and (t<=1)):
	      pre_img=a+t*(b-a)
	      return (pre_img[0],pre_img[1],pre_img[2])
	if (abp_cross[0]!=0):
	   den=(b-a)[2]*p[1]-(b-a)[1]*p[2]
	   num=a[1]*p[2]-a[2]*p[1]
	   t=num/den
	   if ((t>=0) and (t<=1)):
	      pre_img=a+t*(b-a)
	      return (pre_img[0],pre_img[1],pre_img[2])
	print "Error: Attempting to call preim() with a segment and a point not in its perspective projection image."
	return None
    #The function returns the intersection of the ray OP with the support plane of the polygon F.  Thus, the return value is a point in 3D which will project to the point P in the image plane under the perspective projection.
    def supp_F(F,P):
	v_0=np.asarray(F[0])
	v_1=np.asarray(F[1])
	v_2=np.asarray(F[2])
	w_1=v_1-v_0
	w_2=v_2-v_0
	n=np.cross(w_1,w_2)
	p=np.asarray((P[0],P[1],h))	
	den=np.dot(p,n)	
	if (den==0):
	   print "Attempting to call preim_supp_F() with a polygon whose supporting plane is parallel to the ray OP."
	   return None
	num=np.dot(v_0,n)
	t=num/den
	pre_img=(t*P[0],t*P[1],t*h)
	return (pre_img)
    #Given a point in 3D, the function dist(P) returns the distance from the origin to P.
    def dist(P):
       	d=math.sqrt(P[0]**2+P[1]**2+P[2]**2)
       	return (d)
    #Given a list of points in 2D, F_2d, the function diam(F_2d) returns the diameter of the set F_2d, and its minimum and maximum x and y coordinates.
    def diam(F):
	X=[]
	Y=[]
	for f in F:
	    X.append(f[0])
	    Y.append(f[1])
	x_min=min(X)
	y_min=min(Y)
	x_max=max(X)
	y_max=max(Y)
	return (math.sqrt((x_max-x_min)**2+(y_max-y_min)**2),x_min,x_max,y_min,y_max)
    #The function inside(P,F_poly) takes as argumens a point P in 2D and an ordered list of points in 2D which represent the vertices of a polygon, and returns true if P is inside the polygon described by F_poly, and false if P is outside.
    def inside(P,F_poly):
	(D,x_min,x_max,y_min,y_max)=diam(F_poly)
	if (P[0]<x_min or P[0]>x_max): 
	   return False
	if (P[1]<y_min or P[1]>y_max): 
           return False
	secant_ray=(P,(P[0]+D+1,P[1]))
        e_poly=[]						
    	for i in range(len(F_poly)-1):
            e_poly.append((F_poly[i],F_poly[i+1]))
	e_poly.append((F_poly[len(F_poly)-1],F_poly[0]))
	intersection_list=[]
	for i in range(len(e_poly)):
	    int_pt=intersection(secant_ray,e_poly[i])
	    if (int_pt!=None):
	       if (int_pt not in intersection_list):
	          intersection_list.append(intersection(secant_ray,e_poly[i]))
	I=len(intersection_list)%2
	if (I==0): return False #parity count is all that is needed (cf. intersection numbers for polygonal Jordan curves)
	if (I==1): return True
	return
    #F_proj is a list containing the projected vertices of F.
    
    F_proj=proj(F,h)
    #e_proj is a list containing the projections of the edges of F on the image plane as line segments.   
    e_proj=[]						
    for i in range(len(F_proj)-1): e_proj.append((F_proj[i],F_proj[i+1]))
    e_proj.append((F_proj[len(F_proj)-1],F_proj[0]))
    #l_proj is the projection of l to the image plane
    l_proj=proj(l,h)
    #l_e_list is a list containing the endpoints of l_proj and the points of intersection of l_proj and the edges in e_proj
    l_e_list=[]
    l_e_list.append(l_proj[0])
    l_e_list.append(l_proj[1])
    for e in e_proj:
	l_e_inter=intersection(l_proj,e)
        if (l_e_inter!=None): l_e_list.append(l_e_inter)
    delta_y=l_proj[1][1]-l_proj[0][1]
    delta_x=l_proj[1][0]-l_proj[0][0]
    if (abs(delta_x)>=abs(delta_y)):
       def sortKey(coord):
	   return (coord[0])
    else:
       def sortKey(coord):
	   return (coord[1])
    l_e_list_sorted=sorted(l_e_list,key=sortKey)
    #fragments is the list of line subsegments of l to test for occlusion
    fragments=[]
    for i in range(len(l_e_list_sorted)-1):
        fragments.append((l_e_list_sorted[i],l_e_list_sorted[i+1]))
    #fragments_test_pts[] will contain a list containing the midpoints of the line fragments which will be tested for visibility.  If the midpoint is visible, then the fragment will be marked as seen from the optical center.
    
    fragments_test_pts=[]
    for f in fragments:
        mid_pt_f=((f[0][0]+f[1][0])/2,(f[0][1]+f[1][1])/2)
        fragments_test_pts.append(mid_pt_f)
    #clippings will be the list of fragments of l visible from (0,0,0)
    clippings=[]
    for i in range(len(fragments_test_pts)):
        if (inside(fragments_test_pts[i],F_proj)): #the test point is inside the projection of F
           if (dist(preim(l,fragments_test_pts[i]))<dist(supp_F(F,fragments_test_pts[i]))):#the preimage of the test point in l is closer than the preimage of the test point in F (i.e. F doesn't block visibility)
              clippings.append(fragments[i])
        else: clippings.append(fragments[i])
    #The following lines of code will concatenate the line segments which have overlap
    pts=[]

    for i in range(len(clippings)):
        pts.append(clippings[i][0])
        pts.append(clippings[i][1])

    result=[]
    i=0
    while (i<len(pts)):
          a=pts[i]
          j=i+1
          while ((j<len(pts)-2) and (pts[j]==pts[j+1])):
                j=j+2
          b=pts[j]
          i=j+1
          result.append((a,b))
    return result
