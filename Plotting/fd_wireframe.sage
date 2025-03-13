def triangulate(corners):
    base = corners[0]
    triangles = []
    i = 1
    while i + 1 < len(corners):
        triangles.append([base, corners[i], corners[i+1]])
        i += 1
    return triangles

def plot_triangle(hemi, vertices, opacity, frame=False):
    x, y, s, t = var('x y s t')
    g  = lambda x,y: (hemi[1] - (x - hemi[0][0])**2 - (y - hemi[0][1])**2).sqrt()
    x = lambda s,t: vertices[0][0] + s*(vertices[1][0] - vertices[0][0]) + s*t*(vertices[2][0]-vertices[1][0])
    y = lambda s,t: vertices[0][1] + s*(vertices[1][1] - vertices[0][1]) + s*t*(vertices[2][1]-vertices[1][1])
    plot = parametric_plot3d(
        (x(s, t), y(s, t), g(x(s, t), y(s, t))), 
        (0, 1), (0, 1), 
        plot_points=60, 
        color='gray', 
        aspect_ratio=[1,1,1], 
        alpha=opacity,
        zmin = 0,
        dpi=150,
        frame=frame
    )
    
    return plot

def plot_polygon(hemi, corners, opacity=1, frame=False):
    """
    hemi is a tuple with hemi[0] = [x,y] the coordinates of the center and hemi[1] is the square
    radius. Corners is a list of [x,y] coordinates of the projections of the corners of the polyhedron.
    Produces a 3d plot of this single polygon.
    """
    #Verify there are at least 3 corners
    from sage.plot.plot3d.base import Graphics3d
    if len(corners) <= 2:
        return Graphics3d()
    
    #Verify that each corner lies under the hemisphere
    #for pt in corners:
    #    sq_dist_from_center = (hemi[0][0] - pt[0])**2 + (hemi[0][1] - pt[1])**2
    #    print(pt, sq_dist_from_center)
    #    if sq_dist_from_center > hemi[1] + 0.00000001:
    #        raise Exception("Points must lie under hemisphere.")

    #Iterate over triangles that build up the projected polygon
    triangles = triangulate(corners)

    plot = Graphics3d()
    
    for tri in triangles:
        plot += plot_triangle(hemi, tri, opacity, frame)

    return plot

def make_outline_corners(corners):
    outline_corners = []
    outline_corners.append([corners[-1], corners[0]])
    for i in range(len(corners) - 1):
        outline_corners.append([corners[i], corners[i+1]])
    return outline_corners

def plot_arc(hemi, vertices, thickness, frame=False):
    t = var('t')
    g  = lambda x,y: (hemi[1] - (x - hemi[0][0])**2 - (y - hemi[0][1])**2).sqrt()
    x = lambda t: vertices[0][0] + t*(vertices[1][0] - vertices[0][0])
    y = lambda t: vertices[0][1] + t*(vertices[1][1] - vertices[0][1])
    plot = parametric_plot3d(
        (x(t), y(t), g(x(t), y(t))),
        (0,1),
        plot_points=60,
        zmin=0,
        aspect_ratio=[1,1,1],
        thickness=thickness,
        color='black',
        dpi=150,
        frame=frame
    )
    return plot

def plot_outline(hemi, corners, thickness, frame=False):
    outline_corners = make_outline_corners(corners)
    
    from sage.plot.plot3d.base import Graphics3d
    
    plot = Graphics3d()
    for corner_interval in outline_corners:
        plot += plot_arc(hemi, corner_interval, thickness, frame)
        
    return plot

def points_equal_within_error(p1, p2):
    error = 0.000000001
    if (p1[0] - p2[0]).abs() < error and (p1[1] - p2[1]).abs() < error:
        return true
    else:
        return false

def plot_line_segment(p1, p2, thickness, frame=False):
    plot = parametric_plot3d(
        (p1[0] + t*(p2[0] - p1[0]), p1[1] + t*(p2[1] - p1[1]), p1[2] + t*(p2[2] - p1[2])),
        (t,0,1),
        plot_points=60,
        thickness=thickness,
        color='black', 
        dpi=150,
        frame=frame
    )
    return plot

def intersect_3_spheres(hemi0, hemi1, hemi2):
    x0, y0 = hemi0[0]
    r0sq = hemi0[1]
    r0 = sqrt(r0sq)
        
    x1, y1 = hemi1[0]
    r1sq = hemi1[1]
    r1 = sqrt(r1sq)

    #Check if hemi0 and hemi1 intersect
    dist_squared = (x0-x1)**2 + (y0-y1)**2
    #circles intersect if dist_squared <= (r0 + r1)^2 and dist_squared >= (r0 - r1)^2
    if dist_squared > (r0 + r1)**2 or dist_squared < (r0 - r1)**2:
        return []
            
    x2, y2 = hemi2[0]
    r2sq = hemi2[1]
    r2 = sqrt(r2sq)

    #Check if hemi0 and hemi2 intersect
    dist_squared = (x0-x2)**2 + (y0-y2)**2
    if dist_squared > (r0 + r2)**2 or dist_squared < (r0 - r2)**2:
        return []

    #Check if hemi1 and hemi2 intersect
    dist_squared = (x1-x2)**2 + (y1-y2)**2
    if dist_squared > (r1 + r2)**2 or dist_squared < (r1 - r2)**2:
        return []


    #This is the line where hemispheres intersect
    #x*(2*x1-2*x0) + y*(2*y1-2*y0) = (r0sq - r1sq) - (x0**2 - x1**2) - (y0**2 - y1**2)
    #x * x01_coeff + y * y01_coeff = const_term_01
    #x * x02_coeff + y * y02_coeff = const_term_02
    x01_coeff = 2*x1-2*x0
    y01_coeff = 2*y1-2*y0

    x02_coeff = 2*x2-2*x0
    y02_coeff = 2*y2-2*y0

    #If the lines from two intersections are parallel, there is no triple intersection
    are_parallel = (x02_coeff * y01_coeff - y02_coeff * x01_coeff) == 0
    if are_parallel:
        return []

    #Compute the intersection of the two lines
    const_term_01 = (r0sq - r1sq) - (x0**2 - x1**2) - (y0**2 - y1**2)
    const_term_02 = (r0sq - r2sq) - (x0**2 - x2**2) - (y0**2 - y2**2)

    x = 0
    y = 0
    if x01_coeff == 0:
        y = const_term_01 / y01_coeff
        x = (const_term_02 - y * y02_coeff) / x02_coeff
    else:
        #x + y * (y01_coeff/x01_coeff) = (const_term_01/x01_coeff)
        #x = (const_term_01/x01_coeff) - y * (y01_coeff/x01_coeff)
        #((const_term_01/x01_coeff) - y * (y01_coeff/x01_coeff)) * x02_coeff + y * y02_coeff = const_term_02
        y = (const_term_02 - (const_term_01/x01_coeff) * x02_coeff) / (y02_coeff - (y01_coeff/x01_coeff) * x02_coeff)
        x = (const_term_01/x01_coeff) - y * (y01_coeff/x01_coeff)

    #See if the line intersection point actually lies under all 3 hemispheres
    if (x - x0)**2 + (y - y0)**2 <= r0sq and (x - x1)**2 + (y - y1)**2 <= r1sq and (x - x2)**2 + (y - y2)**2 <= r2sq:
        #print("intersect_3_spheres", x, y)
        return [x,y]
    else:
        return []
    
    

def intersect_2_spheres_1_plane(hemi0, hemi1, plane):
    x0, y0 = hemi0[0]
    r0sq = hemi0[1]
    r0 = sqrt(r0sq)
        
    x1, y1 = hemi1[0]
    r1sq = hemi1[1]
    r1 = sqrt(r1sq)

    dist_squared = (x0-x1)**2 + (y0-y1)**2
    #circles intersect if dist_squared <= (r0 + r1)^2 and dist_squared >= (r0 - r1)^2
    if dist_squared > (r0 + r1)**2 or dist_squared < (r0 - r1)**2:
        return []

    x = None
    y = None
    if plane[0] == "y":
        #this is a plane Y = ...
        y = plane[1]
        #means we have to divide by x1-x0
        if x1 == x0:
            return []
        x = (r0sq - r1sq) - (x0**2 - x1**2) - (y0**2 - y1**2) - y*(2*y1 - 2*y0)
        x /= 2*x1 - 2*x0
    elif plane[0] == "x":
        #this is a plane X = ...
        x = plane[1]
        #means we have to divide by y1-y0
        if y1 == y0:
            return []
        y = (r0sq - r1sq) - (x0**2 - x1**2) - (y0**2 - y1**2) - x*(2*x1 - 2*x0)
        y /= 2*y1 - 2*y0
    else:
        raise Exception("Plane argument for \"x or y = c\" should be [x or y, c]")

    if (x - x0)**2 + (y - y0)**2 <= r0sq and (x - x1)**2 + (y - y1)**2 <= r1sq:
        #print("intersect_2_spheres_1_plane", x, y)
        return [x,y]
    else:
        return []

def intersect_1_sphere_2_planes(hemi, plane0, plane1):
    if plane0 == plane1:
        return []
        #raise Exception("Not defined when planes are equal")
        
    x = 0
    y = 0
    if plane0[0] == plane1[0]:
        return []
    elif plane0[0] == "x" and plane1[0] == "y":
        x = plane0[1]
        y = plane1[1]
    else:
        y = plane0[1]
        x = plane1[1]

    if (x - hemi[0][0])**2 + (y - hemi[0][1])**2 <= hemi[1]:
        #print("intersect_1_sphere_2_planes", x,y)
        return [x,y]
    else:
        return []


def plot_scooped_rectangle(plane_axis, plane_value, horiz_min, horiz_max, z_max, hemi, color_function, frame, function_plot_points, colormap):
    u,v = var('u,v')
    X = None
    Y = None
    if plane_axis == "x":
        X = lambda u,v : plane_value
        Y = lambda u,v : u
    elif plane_axis == "y":
        X = lambda u,v : u
        Y = lambda u,v : plane_value
    else:
        raise Exception("plane_axis should be \"x\" or \"y\"")

    hemi_height = lambda u,v : sqrt(hemi[1] - (X(u,v) - hemi[0][0])**2 - (Y(u,v) - hemi[0][1])**2)
    Z = lambda u,v : hemi_height(u,v) + v * (z_max - hemi_height(u,v))

    cf = lambda u,v : color_function(X(u,v), Y(u,v), Z(u,v))
    plot = parametric_plot3d(
        [X, Y, Z],
        (u,horiz_min,horiz_max),
        (v,0,1),
        plot_points=function_plot_points,
        color = (cf,colormap),
        frame=frame
        )
    return plot
	
#This is the "does it all" function.
#Parameter functionality is explained in the examples
#Always evaluate this cell

def fd_wireframe(d, verbose=0, opacity=1, thickness=2,
                 plot_function=False,
                 f=0, plane_axis="x", plane_value="0",
                 animate=False,
                 values_list=[],
                 zmax=3,
                 frame=False,
                 function_plot_points=200,
                 function_colormap = None,
                 #browser="/usr/bin/firefox"):
                 browser="/opt/google/chrome/chrome"):
    """
    Returns a 3d plot of a fundamental domain for the Bianchi group
    for x^2 + d
    """

    from utils import (make_M_alphas,
                       make_poly_from_edges,
                       poly_gl2_orbit_reps,
                       aas_triangle_gl2_orbit_reps,
                       oriented_faces, face_index,
                       polyhedron_relation,
                       polygon_parameters,
                       std_poly,
                       is_poly_principal)
    from H3 import (make_k, 
                    alphas_sigmas_from_file,
                    output_alphas_sigmas,
                    all_polyhedra)

    if plot_function == True:
        if plane_axis != "x" and plane_axis != "y":
            raise Exception("plane_axis should equal \"x\" or \"y\"")
            
            
    colormap = colormaps.viridis
    if not function_colormap is None:
        colormap = function_colormap
    
    kdata = make_k(d)
    k = kdata['k']
    if verbose:
        print("Field: {}".format(k))
        print("Discriminant: {}".format(k.discriminant()))
        print("Class number: {}".format(k.class_number()))

    # NB we may have precomputed alphas and sigmas and put the data
    # into a geodata file but not yet computed the tessellation
    geodata_file = f"geodata_{d}.dat"
    try:
        alphas, sigmas, plus_pairs, minus_pairs, fours = alphas_sigmas_from_file(kdata)
        if verbose:
            print("using precomputed alphas")
    except ValueError:
        if verbose:
            print("no precomputed alphas, computing from scratch")
        alphas, sigmas, plus_pairs, minus_pairs, fours = alpha_sigma_data(kdata, verbose=verbose)
    with open(geodata_file, 'w') as geout:
        output_alphas_sigmas(kdata, plus_pairs, minus_pairs, fours, sigmas, geout)

    M_alphas, alpha_inv = make_M_alphas(alphas)
    if verbose:
        print("{} alphas".format(len(alphas)))
        print("{} sigmas".format(len(sigmas)))

    polyhedra, hemis = all_polyhedra(k, alphas, verbose>1)
    npoly = len(polyhedra)
    poly = "polyhedra" if npoly>1 else "polyhedron"
    print(f"Tessellation has {npoly} {poly}:")
    print(polyhedra)

    print("plotting fundamental domain wireframe")
    from sage.misc.viewer import viewer
    viewer.browser(browser)

    from sage.plot.plot3d.base import Graphics3d
    
    hemispheres = []
    for H in hemis:
        x, y = kdata['emb'](H[0])
        rsq = H[1]
        r = sqrt(rsq)
        hemispheres.append([[x,y], rsq])

    xmax = 1/2
    xmin = -xmax
    ymax = kdata['Ymax']
    ymin = -ymax

    if plot_function == True and animate == False:
        if plane_axis == "x":
            if plane_value < xmin or xmax < plane_value:
                raise Exception("plane_value should be between min and max values")
        if plane_axis == "y":
            if plane_value < ymin or ymax < plane_value:
                raise Exception("plane_value should be between min and max values")
    elif plot_function == True and animate == True:
        if plane_axis == "x":
            for plane_value in values_list:
                if plane_value < xmin or xmax < plane_value:
                    raise Exception("plane_value should be between min and max values")
        if plane_axis == "y":
            for plane_value in values_list:
                if plane_value < ymin or ymax < plane_value:
                    raise Exception("plane_value should be between min and max values")
    
    boundary_planes = []
    boundary_planes.append(["x",xmax])
    boundary_planes.append(["x",xmin])
    boundary_planes.append(["y",ymax])
    boundary_planes.append(["y",ymin])

    plane_index = len(hemispheres)

    surfaces = hemispheres + boundary_planes
    
    hemi_corners_list = []
    for i in range(len(hemispheres)):
        hemi_corners_list.append([])
        
        for j in range(len(surfaces)):
            for k in range(j + 1, len(surfaces)):
                corner = []
                if i == j or i == k or j == k:
                    continue
                elif j >= plane_index and k >= plane_index:
                    #i is sphere, j is plane, k is plane
                    corner = intersect_1_sphere_2_planes(surfaces[i], surfaces[j], surfaces[k])
                    
                elif j >= plane_index and k < plane_index:
                    #i is sphere, j is plane, k is sphere
                    corner = intersect_2_spheres_1_plane(surfaces[i], surfaces[k], surfaces[j])
                        
                elif j < plane_index and k >= plane_index:
                    #i is sphere, j is sphere, k is plane
                    corner = intersect_2_spheres_1_plane(surfaces[i], surfaces[j], surfaces[k])
                        
                else:
                    #i is sphere, j is sphere, k is sphere
                    corner = intersect_3_spheres(hemispheres[i], hemispheres[j], hemispheres[k])
                
                if corner:
                     if xmin <= corner[0] and corner[0] <= xmax and ymin <= corner[1] and corner[1] <= ymax:
                        hemi_corners_list[i].append(corner)

    #Clean up list:
    #Convert to floating point
    for i in range(len(hemi_corners_list)):
        for j in range(len(hemi_corners_list[i])):
            hemi_corners_list[i][j] = [n(hemi_corners_list[i][j][0]), n(hemi_corners_list[i][j][1])]

    #Remove duplicates that are within floating point error
    for i in range(len(hemi_corners_list)):
        temp = []
        for j in range(len(hemi_corners_list[i])):
            pt = hemi_corners_list[i][j]
            #if pt is not in temp, add it
            in_temp = False
            for k in range(len(temp)):
                pt2 = temp[k]
                if points_equal_within_error(pt, pt2):
                    in_temp = True
                    break
            if not in_temp:
                temp.append(pt)
            
        hemi_corners_list[i] = temp
    
    #Delete corners which are strictly underneath another hemisphere
    #I expect JC will figure out a way so that I can get rid of this step
    for i in range(len(hemi_corners_list)):
        temp = []
        this_hemi = hemispheres[i]
        for corner in hemi_corners_list[i]:
            zsq = (this_hemi[1] - (corner[0] - this_hemi[0][0])**2 - (corner[1] - this_hemi[0][1])**2)
            remove_this_corner = False
            for j in range(len(hemi_corners_list)):
                if i == j:
                    continue
                hemi = hemispheres[j]

                if (corner[0] - hemi[0][0])**2 + (corner[1] - hemi[0][1])**2 + zsq < hemi[1] - 0.000001:
                    remove_this_corner = True
                    break
            if not remove_this_corner:
                temp.append(corner)
        hemi_corners_list[i] = temp


    #Sort each list by angle with respect to their centroid
    for i in range(len(hemi_corners_list)):
        x_average = 0
        y_average = 0
        for j in range(len(hemi_corners_list[i])):
            x_average += hemi_corners_list[i][j][0]
            y_average += hemi_corners_list[i][j][1]
        x_average /= len(hemi_corners_list[i])
        y_average /= len(hemi_corners_list[i])
        
        hemi_corners_list[i].sort(key=lambda pt: atan2(pt[1] - y_average, pt[0] - x_average)) 

    #Begin plotting
    plot = Graphics3d()
    if opacity > 0:
        #Plot polygons
        for i in range(len(hemispheres)):
            plot += plot_polygon(hemispheres[i], hemi_corners_list[i], opacity, frame=frame)
        

    #Plot polygon outlines
    for i in range(len(hemispheres)):
        plot += plot_outline(hemispheres[i], hemi_corners_list[i], thickness, frame=frame)

    #Plot non-curved lines of fundamental domain
    #makes the assumption that the vertical lines will have the same lowest z-coordinate because it's an FD!
    corner_height = 0
    for hemi in hemispheres:
        if (xmax - hemi[0][0])**2 + (ymax - hemi[0][1])**2 <= hemi[1]:
            corner_height = (hemi[1] - (xmax - hemi[0][0])**2 - (ymax - hemi[0][1])**2).sqrt()
            break
    
    top_height = zmax
    plot += plot_line_segment((xmin, ymin, corner_height), (xmin, ymin, top_height), thickness=thickness, frame=frame)
    plot += plot_line_segment((xmin, ymax, corner_height), (xmin, ymax, top_height), thickness=thickness, frame=frame)
    plot += plot_line_segment((xmax, ymin, corner_height), (xmax, ymin, top_height), thickness=thickness, frame=frame)
    plot += plot_line_segment((xmax, ymax, corner_height), (xmax, ymax, top_height), thickness=thickness, frame=frame)
    plot += plot_line_segment((xmin, ymin, top_height), (xmax, ymin, top_height), thickness=thickness, frame=frame)
    plot += plot_line_segment((xmax, ymin, top_height), (xmax, ymax, top_height), thickness=thickness, frame=frame)
    plot += plot_line_segment((xmax, ymax, top_height), (xmin, ymax, top_height), thickness=thickness, frame=frame)
    plot += plot_line_segment((xmin, ymax, top_height), (xmin, ymin, top_height), thickness=thickness, frame=frame)

    if plot_function == False:
        return plot

    #Compute where this plane intersects hemispheres
    plots = []

    if animate == False:
        c = plane_value
        #Set plane to represent the plane that we're plotting on
        plane = [plane_axis, c]
    
        #restrict to only the hemispheres that CAN intersect with the plane
        restricted_hemispheres = []
        for hemi in hemispheres:
            if plane_axis == "x":
                if (c - hemi[0][0])**2 <= hemi[1]:
                    restricted_hemispheres.append(hemi)
            else:
                if (c - hemi[0][1])**2 <= hemi[1]:
                    restricted_hemispheres.append(hemi)
    
        #compute intersections with the plane and two hemispheres
        hemispheres_and_intervals = []
        for i in range(len(restricted_hemispheres)):
            hemi0 = restricted_hemispheres[i]
            other_coords = []
            for j in range(len(restricted_hemispheres)):
                if i == j:
                    continue
                hemi1 = restricted_hemispheres[j]
                xy_coords = intersect_2_spheres_1_plane(hemi0, hemi1, plane)
                    
                if xy_coords:
                    if plane_axis == "x":
                        other_coords.append(xy_coords[1])
                    else:
                        other_coords.append(xy_coords[0])
                        
            hemispheres_and_intervals.append([hemi0, other_coords])
    
        #compute intersections with the plane, a hemisphere, and another plane
        for i in range(len(restricted_hemispheres)):
            hemi = restricted_hemispheres[i]
            for boundary_plane in boundary_planes:
                xy_coords = intersect_1_sphere_2_planes(hemi, plane, boundary_plane)
    
                if xy_coords:
                    if plane_axis == "x":
                        hemispheres_and_intervals[i][1].append(xy_coords[1])
                    else:
                        hemispheres_and_intervals[i][1].append(xy_coords[0])
    
        #delete extraneous hemispheres
        extraneous_hai_deleted = []
        for i in range(len(hemispheres_and_intervals)):
            hai0 = hemispheres_and_intervals[i]
            hemi0 = hai0[0]
            new_hai = [hemi0,[]]
            
            for other_coord in hai0[1]:
                keep_this_coord = True
                for j in range(len(hemispheres_and_intervals)):
                    if i == j:
                        continue
                    hai1 = hemispheres_and_intervals[j]
                    hemi1 = hai1[0]
        
                    if plane_axis == "x":
                        if (c - hemi1[0][0])**2 + (other_coord - hemi1[0][1])**2 - (c - hemi0[0][0])**2 - (other_coord - hemi0[0][1])**2  < hemi1[1] - hemi0[1] - 0.000001:
                            keep_this_coord = False
                            break
                    else:
                        if (other_coord - hemi1[0][0])**2 + (c - hemi1[0][1])**2 - (other_coord - hemi0[0][0])**2 - (c - hemi0[0][1])**2  < hemi1[1] - hemi0[1] - 0.000001:
                            keep_this_coord = False
                            break
                if keep_this_coord:
                    new_hai[1].append(other_coord)
            extraneous_hai_deleted.append(new_hai)
            
        final_hemis_and_values = []
        for i in range(len(extraneous_hai_deleted)):
            hai = extraneous_hai_deleted[i]
            if len(hai[1]) >= 2:
                if min(hai[1]) != max(hai[1]):
                    if plane_axis == "x":
                        left = max(ymin, min(hai[1]))
                        right = min(ymax, max(hai[1]))
                        final_hemis_and_values.append([hai[0], [left, right]])
                    else:
                        left = max(xmin, min(hai[1]))
                        right = min(xmax, max(hai[1]))
                        final_hemis_and_values.append([hai[0], [left, right]])

        plot_this_plane = plot
        for [hemi, interval] in final_hemis_and_values:
            plot_this_plane += plot_scooped_rectangle(plane_axis=plane_axis, 
                                           plane_value=c, 
                                           horiz_min=interval[0], 
                                           horiz_max=interval[1], 
                                           z_max=top_height, 
                                           hemi=hemi, 
                                           color_function=f,
                                           frame=frame,
                                           function_plot_points=function_plot_points,
                                           colormap=colormap)

        #Options: https://doc.sagemath.org/html/en/reference/plot3d/sage/plot/plot3d/base.html#sage.plot.plot3d.base.Graphics3d.show
        cam_pos = (2,2,3)
        view_target = (0,0,1)
        view_dir = (cam_pos[0] - view_target[0], cam_pos[1] - view_target[1], cam_pos[2] - view_target[2])
        import time
        ts = time.time()
        import datetime
        fname = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S') 
        plot_this_plane.save(fname + ".png", 
                     viewer='tachyon',
                     shade='medium',
                     zoom = 1.5,
                     figsize=[10,10],
                     aspect_ratio=[1,1,1], 
                     antialiasing=8,
                    look_at=(0,0,1))
        return plot_this_plane
    else:
        import time
        ts = time.time()
        import datetime
        fname = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        import os
        folder = "frames_" + str(d) + "_" + fname
        try:
            os.mkdir(folder)
        except Exception as e:
            print(f"An error occurred: {e}")
        
        
        count = 0
        for c in values_list:
            #Set plane to represent the plane that we're plotting on
            plane = [plane_axis, c]
        
            #restrict to only the hemispheres that CAN intersect with the plane
            restricted_hemispheres = []
            for hemi in hemispheres:
                if plane_axis == "x":
                    if (c - hemi[0][0])**2 <= hemi[1]:
                        restricted_hemispheres.append(hemi)
                else:
                    if (c - hemi[0][1])**2 <= hemi[1]:
                        restricted_hemispheres.append(hemi)
        
            #compute intersections with the plane and two hemispheres
            hemispheres_and_intervals = []
            for i in range(len(restricted_hemispheres)):
                hemi0 = restricted_hemispheres[i]
                other_coords = []
                for j in range(len(restricted_hemispheres)):
                    if i == j:
                        continue
                    hemi1 = restricted_hemispheres[j]
                    xy_coords = intersect_2_spheres_1_plane(hemi0, hemi1, plane)
                        
                    if xy_coords:
                        if plane_axis == "x":
                            other_coords.append(xy_coords[1])
                        else:
                            other_coords.append(xy_coords[0])
                            
                hemispheres_and_intervals.append([hemi0, other_coords])
        
            #compute intersections with the plane, a hemisphere, and another plane
            for i in range(len(restricted_hemispheres)):
                hemi = restricted_hemispheres[i]
                for boundary_plane in boundary_planes:
                    xy_coords = intersect_1_sphere_2_planes(hemi, plane, boundary_plane)
        
                    if xy_coords:
                        if plane_axis == "x":
                            hemispheres_and_intervals[i][1].append(xy_coords[1])
                        else:
                            hemispheres_and_intervals[i][1].append(xy_coords[0])
        
            #delete extraneous hemispheres
            extraneous_hai_deleted = []
            for i in range(len(hemispheres_and_intervals)):
                hai0 = hemispheres_and_intervals[i]
                hemi0 = hai0[0]
                new_hai = [hemi0,[]]
                
                for other_coord in hai0[1]:
                    keep_this_coord = True
                    for j in range(len(hemispheres_and_intervals)):
                        if i == j:
                            continue
                        hai1 = hemispheres_and_intervals[j]
                        hemi1 = hai1[0]
            
                        if plane_axis == "x":
                            if (c - hemi1[0][0])**2 + (other_coord - hemi1[0][1])**2 - (c - hemi0[0][0])**2 - (other_coord - hemi0[0][1])**2  < hemi1[1] - hemi0[1] - 0.000001:
                                keep_this_coord = False
                                break
                        else:
                            if (other_coord - hemi1[0][0])**2 + (c - hemi1[0][1])**2 - (other_coord - hemi0[0][0])**2 - (c - hemi0[0][1])**2  < hemi1[1] - hemi0[1] - 0.000001:
                                keep_this_coord = False
                                break
                    if keep_this_coord:
                        new_hai[1].append(other_coord)
                extraneous_hai_deleted.append(new_hai)
                
            final_hemis_and_values = []
            for i in range(len(extraneous_hai_deleted)):
                hai = extraneous_hai_deleted[i]
                if len(hai[1]) >= 2:
                    if min(hai[1]) != max(hai[1]):
                        if plane_axis == "x":
                            left = max(ymin, min(hai[1]))
                            right = min(ymax, max(hai[1]))
                            final_hemis_and_values.append([hai[0], [left, right]])
                        else:
                            left = max(xmin, min(hai[1]))
                            right = min(xmax, max(hai[1]))
                            final_hemis_and_values.append([hai[0], [left, right]])
    
            plot_this_plane = plot
            for [hemi, interval] in final_hemis_and_values:
                plot_this_plane += plot_scooped_rectangle(plane_axis=plane_axis, 
                                               plane_value=c, 
                                               horiz_min=interval[0], 
                                               horiz_max=interval[1], 
                                               z_max=top_height, 
                                               hemi=hemi, 
                                               color_function=f,
                                               frame=frame,
                                               function_plot_points=function_plot_points,
                                               colormap=colormap)
    
            #Options: https://doc.sagemath.org/html/en/reference/plot3d/sage/plot/plot3d/base.html#sage.plot.plot3d.base.Graphics3d.show
            cam_pos = (2,2,3)
            view_target = (0,0,1)
            view_dir = (cam_pos[0] - view_target[0], cam_pos[1] - view_target[1], cam_pos[2] - view_target[2]) 
            file_number_str = str(count)
            while len(file_number_str) < 3:
                file_number_str = "0" + file_number_str

            plot_this_plane.save(folder + "/" + file_number_str + ".png", 
                         viewer='tachyon',
                         shade='medium',
                         zoom = 1.5,
                         figsize=[10,10],
                         aspect_ratio=[1,1,1], 
                         antialiasing=8,
                        look_at=(0,0,1))
            print("frame " + str(count) + " done")
            count += 1
        return plot_this_plane        