/*
-----------------------------------------------------------------------------
This source file has been developed within the scope of the
Technical Director course at Filmakademie Baden-Wuerttemberg.
http://technicaldirector.de
    
Written by Juraj Tomori.
Copyright (c) 2019 Animationsinstitut of Filmakademie Baden-Wuerttemberg
-----------------------------------------------------------------------------
*/

#ifndef VFT_SUBD
#define VFT_SUBD

// a function for constructing new primitives, to be used in subdivision functions
function void construct_primitives(int prim_pts_idx[], new_edge_pts[]; int new_face_pt, ordering_mode, primnum)
{
    removeprim(0, primnum, 0);

    if (ordering_mode == 0)
    {
        // uniform
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", new_edge_pts[0], prim_pts_idx[1], new_edge_pts[1], new_face_pt);
        addprim(0, "poly", new_face_pt, new_edge_pts[1], prim_pts_idx[2], new_edge_pts[2]);
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[2], prim_pts_idx[3]);
    }
    else if (ordering_mode == 1)
    {
        // mirror x y
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", prim_pts_idx[1], new_edge_pts[0], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", prim_pts_idx[2], new_edge_pts[2], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", prim_pts_idx[3], new_edge_pts[2], new_face_pt, new_edge_pts[3]);
    }
    else if (ordering_mode == 2)
    {
        // mirror y
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", prim_pts_idx[1], new_edge_pts[0], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", new_edge_pts[1], new_face_pt, new_edge_pts[2], prim_pts_idx[2]);
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[2], prim_pts_idx[3]);
    }
    else if (ordering_mode == 3)
    {
        // mirror x
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", new_edge_pts[0], prim_pts_idx[1], new_edge_pts[1], new_face_pt);
        addprim(0, "poly", new_edge_pts[2], prim_pts_idx[2], new_edge_pts[1], new_face_pt);
        addprim(0, "poly", prim_pts_idx[3], new_edge_pts[2], new_face_pt, new_edge_pts[3]);
    }
    else if (ordering_mode == 4)
    {
        // cyclic
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", prim_pts_idx[1], new_edge_pts[1], new_face_pt, new_edge_pts[0]);
        addprim(0, "poly", prim_pts_idx[2], new_edge_pts[2], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", prim_pts_idx[3], new_edge_pts[3], new_face_pt, new_edge_pts[2]);
    }
    else if (ordering_mode == 5)
    {
        // opposite
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[0], prim_pts_idx[0]);
        addprim(0, "poly", new_edge_pts[1], new_face_pt, new_edge_pts[0], prim_pts_idx[1]);
        addprim(0, "poly", new_edge_pts[1], new_face_pt, new_edge_pts[2], prim_pts_idx[2]);
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[2], prim_pts_idx[3]);
    }
}

// returns point number of a vertex in a primitive
// primnum - primitive containing a vertex
// vertexnum - vertex number (relative to the primitive, starting at 0)
function int vertexprimpoint(int geo, primnum, vertexnum)
{
    return vertexpoint(geo, vertexindex(geo, primnum, vertexnum) );
}

// a struct containing weights and parameters - to be used for transferring all the weights into functions
struct quad_face_split_parms
{
    float face_weights[] = {1.0, 1.0, 1.0, 1.0};
    float face_weight_offset = 0.0;
    float edge_weights[] = {1.0, 1.0};
    float edge_weight_offset = 0.0;
    int ordering_mode = 0;
}

// basic quad subdivision
function void quad_face_split(int primnum; quad_face_split_parms weights)
{
    ////// face point
    
    // find pints in a primitive
    int prim_pts_idx[] = primpoints(0, primnum);
    vector new_face_P = {0, 0, 0};

    // interpolate attribute
    for (int i = 0; i < len(prim_pts_idx); i++)
    {
        new_face_P += vector(point(0, "P", prim_pts_idx[i])) * weights.face_weights[i];
    }

    new_face_P /= sum(weights.face_weights) + weights.face_weight_offset;

    // add new point
    int new_face_pt = addpoint(0, new_face_P);

    ////// edge points

    vector2 edge_vertices[] = array( {0, 1}, {1, 2}, {2, 3}, {3, 0} );
    int new_edge_pts[];

    foreach (vector2 edge; edge_vertices)
    {
        // get point numbers of vertices in the edge
        int pt_1 = vertexprimpoint(0, primnum, (int)edge[0]);
        int pt_2 = vertexprimpoint(0, primnum, (int)edge[1]);

        // interpolate attribute
        vector pt_1_P = point(0, "P", pt_1);
        vector pt_2_P = point(0, "P", pt_2);

        vector new_edge_P = pt_1_P * weights.edge_weights[0] + pt_2_P * weights.edge_weights[1];
        new_edge_P /= sum(weights.edge_weights) + weights.edge_weight_offset;

        // add new point
        append(new_edge_pts, addpoint(0, new_edge_P));
    }

    ////// construct primitives
    construct_primitives(prim_pts_idx, new_edge_pts, new_face_pt, weights.ordering_mode, primnum);
}

struct catmull_clark_quad_parms
{
    float new_face_pt_weights[] = {1.0, 1.0, 1.0, 1.0};
    float new_face_pt_weights_offset = 0.0;
    float new_edge_pt_weights[] = {6.0, 6.0, 1.0, 1.0, 1.0, 1.0};
    float new_edge_pt_weights_offset = 0.0;
    float pt_3_pt_weights[] = {15.0, 6.0, 6.0, 6.0, 1.0, 1.0, 1.0};
    float pt_3_pt_weights_offset = 0.0;
    float pt_4_pt_weights[] = {36.0, 6.0, 6.0, 6.0, 6.0, 1.0, 1.0, 1.0, 1.0};
    float pt_4_pt_weights_offset = 0.0;
    float pt_5_pt_weights[] = {65.0, 6.0, 6.0, 6.0, 6.0, 6.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float pt_5_pt_weights_offset = 0.0;
    
    int ordering_mode = 0;
}

// catmull clark subdivision on quad meshes
function void catmull_clark_quad(int primnum; catmull_clark_quad_parms weights)
{
    ////// face point
    
    // find pints in a primitive
    int prim_pts_idx[] = primpoints(0, primnum);
    vector new_face_P = {0, 0, 0};

    // interpolate attribute
    for (int i = 0; i < len(prim_pts_idx); i++)
    {
        new_face_P += vector(point(0, "P", prim_pts_idx[i])) * weights.new_face_pt_weights[i];
    }

    new_face_P /= sum(weights.new_face_pt_weights) + weights.new_face_pt_weights_offset;

    // add new point
    int new_face_pt = addpoint(0, new_face_P);

    ////// edge points
    
    vector2 edge_vertices[] = array( {0, 1}, {1, 2}, {2, 3}, {3, 0} );
    int new_edge_pts[];

    foreach (vector2 edge; edge_vertices)
    {
        // get point numbers of vertices in the edge
        int pt_1 = vertexprimpoint(0, primnum, (int)edge[0]);
        int pt_2 = vertexprimpoint(0, primnum, (int)edge[1]);

        // find all primitives that contain the two points
        int prims_1[] = pointprims(0, pt_1);
        int prims_2[] = pointprims(0, pt_2);
        int prims_shared[];

        // find shared primitives
        foreach (int pt_prim; prims_1)
        {
            int found_idx = find(prims_2, pt_prim);
            if (found_idx >= 0)
                append(prims_shared, pt_prim);
        }

        if (len(prims_shared) > 2)
            warning("Edge | Found more than two shared primitives");

        // if edge is shared between two primitives
        if (len(prims_shared) == 2)
        {
            int edge_prims_pts[];
            append(edge_prims_pts, primpoints(0, prims_shared[0]));
            append(edge_prims_pts, primpoints(0, prims_shared[1]));

            removevalue(edge_prims_pts, pt_1);
            removevalue(edge_prims_pts, pt_1);
            removevalue(edge_prims_pts, pt_2);
            removevalue(edge_prims_pts, pt_2);
            edge_prims_pts = sort(edge_prims_pts);

            if (len(edge_prims_pts) != 4)
                warning("Edge | There are not 4 unshared points between two primitives");
            
            int influence_pt_nums[];
            append(influence_pt_nums, sort(array(pt_1, pt_2)));
            append(influence_pt_nums, edge_prims_pts);

            if (len(influence_pt_nums) != 6)
                warning("Edge | Length of array with influencing points is not 6");

            vector new_edge_P = {0, 0, 0};

            // interpolate attribute
            for (int i = 0; i < 6; i++)
            {
                new_edge_P += vector(point(0, "P", influence_pt_nums[i])) * weights.new_edge_pt_weights[i];
            }

            new_edge_P /= sum(weights.new_edge_pt_weights) + weights.new_edge_pt_weights_offset;

            // add new point
            append(new_edge_pts, addpoint(0, new_edge_P));
        }

        // if edge is not shared between primitives
        else if (len(prims_shared) == 1)
        {
            // interpolate attribute
            vector pt_1_P = point(0, "P", pt_1);
            vector pt_2_P = point(0, "P", pt_2);

            vector new_edge_P = pt_1_P * weights.new_edge_pt_weights[0] + pt_2_P * weights.new_edge_pt_weights[1];
            new_edge_P /= weights.new_edge_pt_weights[0] + weights.new_edge_pt_weights[1] + weights.new_edge_pt_weights_offset;

            // add new point
            append(new_edge_pts, addpoint(0, new_edge_P));
        }
    }

    ////// construct primitives

    construct_primitives(prim_pts_idx, new_edge_pts, new_face_pt, weights.ordering_mode, primnum);

    ////// original vertices

    foreach (int prim_pt; prim_pts_idx)
    {
        int prim_pt_prims[] = sort(pointprims(0, prim_pt));

        // only distort vertex if the current primitive has the highest number of primitives that belong to the vertex - to avoid duplicate translations
        if (primnum == prim_pt_prims[-1])
        {
            int influence_pt_nums[];
            append(influence_pt_nums, prim_pt); // myself

            int prim_pt_neigbors[] = sort(neighbours(0, prim_pt));
            append(influence_pt_nums, prim_pt_neigbors); // neighbors
            int prim_pt_valence = len(prim_pt_neigbors);

            int prim_pt_distant_neighbors[];

            foreach(int prim_pt_prim; prim_pt_prims) // for each primitive that belongs to iterated point
            {
                int prim_pt_prim_pts[] = primpoints(0, prim_pt_prim);

                removevalue(prim_pt_prim_pts, prim_pt); // remove itself

                foreach (int prim_pt_neigbor; prim_pt_neigbors)
                {
                    removevalue(prim_pt_prim_pts, prim_pt_neigbor); // remove neighbors
                }

                if (len(prim_pt_prim_pts) != 1)
                    warning("Original vertex | After removing center and neighbor vertices, there is not exactly one point left");

                append(influence_pt_nums, prim_pt_prim_pts[0]); // distant neighbors
            }

            vector new_orig_vertex_P = {0, 0, 0};

            //if (prim_pt_valence > 5)
            //        error("Original vertex | Subdivision is supported only for geometry with points of maximum valence 5");
            
            // interpolate attribute
            int neighbor_counter = 0;
            foreach (int influence_pt_num; influence_pt_nums)
            {
                if (prim_pt_valence == 3)
                {
                    new_orig_vertex_P += vector(point(0, "P", influence_pt_num)) * weights.pt_3_pt_weights[neighbor_counter];
                }
                else if (prim_pt_valence == 4)
                {
                    new_orig_vertex_P += vector(point(0, "P", influence_pt_num)) * weights.pt_4_pt_weights[neighbor_counter];
                }
                else if (prim_pt_valence == 5)
                {
                    new_orig_vertex_P += vector(point(0, "P", influence_pt_num)) * weights.pt_5_pt_weights[neighbor_counter];
                }

                neighbor_counter++;
            }

            if (prim_pt_valence == 3)
                new_orig_vertex_P /= sum(weights.pt_3_pt_weights) + weights.pt_3_pt_weights_offset;
            else if (prim_pt_valence == 4)
                new_orig_vertex_P /= sum(weights.pt_4_pt_weights) + weights.pt_4_pt_weights_offset;
            else if (prim_pt_valence == 5)
                new_orig_vertex_P /= sum(weights.pt_5_pt_weights) + weights.pt_5_pt_weights_offset;
            
            // translate
            if (prim_pt_valence >= 3)
                setpointattrib(0, "P", prim_pt, new_orig_vertex_P, "set");
        }
    }
}

#endif
