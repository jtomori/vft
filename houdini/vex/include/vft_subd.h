#ifndef VFT_SUBD
#define VFT_SUBD

// a struct containing weights and parameters - to be used for transferring all the weights into functions
struct quad_face_split_parms
{
    float face_weights[] = {1.0, 1.0, 1.0, 1.0};
    float face_weight_offset = 0.0;
    float edge_weights[] = {1.0, 1.0};
    float edge_weight_offset = 0.0;
    int ordering_mode = 0;
}

// returns point number of a vertex in a primitive
// primnum - primitive containing a vertex
// vertexnum - vertex number (relative to the primitive, starting at 0)
function int vertexprimpoint(int geo, primnum, vertexnum)
{
    return vertexpoint(geo, vertexindex(geo, primnum, vertexnum) );
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

    removeprim(0, primnum, 0);

    if (weights.ordering_mode == 0)
    {
        // uniform
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", new_edge_pts[0], prim_pts_idx[1], new_edge_pts[1], new_face_pt);
        addprim(0, "poly", new_face_pt, new_edge_pts[1], prim_pts_idx[2], new_edge_pts[2]);
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[2], prim_pts_idx[3]);
    }
    else if (weights.ordering_mode == 1)
    {
        // mirror x y
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", prim_pts_idx[1], new_edge_pts[0], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", prim_pts_idx[2], new_edge_pts[2], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", prim_pts_idx[3], new_edge_pts[2], new_face_pt, new_edge_pts[3]);
    }
    else if (weights.ordering_mode == 2)
    {
        // mirror y
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", prim_pts_idx[1], new_edge_pts[0], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", new_edge_pts[1], new_face_pt, new_edge_pts[2], prim_pts_idx[2]);
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[2], prim_pts_idx[3]);
    }
    else if (weights.ordering_mode == 3)
    {
        // mirror x
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", new_edge_pts[0], prim_pts_idx[1], new_edge_pts[1], new_face_pt);
        addprim(0, "poly", new_edge_pts[2], prim_pts_idx[2], new_edge_pts[1], new_face_pt);
        addprim(0, "poly", prim_pts_idx[3], new_edge_pts[2], new_face_pt, new_edge_pts[3]);
    }
    else if (weights.ordering_mode == 4)
    {
        // cyclic
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", prim_pts_idx[1], new_edge_pts[1], new_face_pt, new_edge_pts[0]);
        addprim(0, "poly", prim_pts_idx[2], new_edge_pts[2], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", prim_pts_idx[3], new_edge_pts[3], new_face_pt, new_edge_pts[2]);
    }
    else if (weights.ordering_mode == 5)
    {
        // opposite
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[0], prim_pts_idx[0]);
        addprim(0, "poly", new_edge_pts[1], new_face_pt, new_edge_pts[0], prim_pts_idx[1]);
        addprim(0, "poly", new_edge_pts[1], new_face_pt, new_edge_pts[2], prim_pts_idx[2]);
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[2], prim_pts_idx[3]);
    }
}

struct catmull_clark_quad_parms
{
    float new_face_pt_weights[] = {1.0, 1.0, 1.0, 1.0};
    float new_edge_pt_weights[] = {6.0, 6.0, 1.0, 1.0, 1.0, 1.0};
    float pt_3_pt_weights[] = {15.0, 6.0, 6.0, 6.0, 1.0, 1.0, 1.0};
    float pt_4_pt_weights[] = {36.0, 6.0, 6.0, 6.0, 6.0, 1.0, 1.0, 1.0, 1.0};
    float pt_5_pt_weights[] = {65.0, 6.0, 6.0, 6.0, 6.0, 6.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    
    int ordering_mode = 0;
}

// catmull clark subdivision on quad meshes
function void catmull_clark_quad(int primnum; catmull_clark_quad_parms weights)
{
    ////// face point
    
    // find pints in a primitive
    int prim_pts[] = primpoints(0, primnum);
    vector new_face_P = {0, 0, 0};

    // interpolate attribute
    for (int i = 0; i < len(prim_pts); i++)
    {
        new_face_P += vector(point(0, "P", prim_pts[i])) * weights.new_face_pt_weights[i];
    }

    new_face_P /= sum(weights.new_face_pt_weights);

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

        // find common neighbors
        int nb_pts[];
        append(nb_pts, neighbours(0, pt_1));
        append(nb_pts, neighbours(0, pt_2));

        // remove source primitive points
        foreach (int pt; prim_pts)
        {
            removevalue(nb_pts, pt);
        }

        // find duplicates / shared neighbors
        sort(nb_pts);
        int nb_prev = nb_pts[0];
        for (int i = 1; i < len(nb_pts); i++)
        {

        }

        /*
        // interpolate attribute
        vector pt_1_P = point(0, "P", pt_1);
        vector pt_2_P = point(0, "P", pt_2);

        vector new_edge_P = pt_1_P * weights.edge_weights[0] + pt_2_P * weights.edge_weights[1];
        new_edge_P /= sum(weights.edge_weights);

        // add new point
        append(new_edge_pts, addpoint(0, new_edge_P));
        */
    }
    

    ////// construct primitives
    /*
    removeprim(0, primnum, 0);

    if (weights.ordering_mode == 0)
    {
        // uniform
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", new_edge_pts[0], prim_pts_idx[1], new_edge_pts[1], new_face_pt);
        addprim(0, "poly", new_face_pt, new_edge_pts[1], prim_pts_idx[2], new_edge_pts[2]);
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[2], prim_pts_idx[3]);
    }
    else if (weights.ordering_mode == 1)
    {
        // mirror x y
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", prim_pts_idx[1], new_edge_pts[0], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", prim_pts_idx[2], new_edge_pts[2], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", prim_pts_idx[3], new_edge_pts[2], new_face_pt, new_edge_pts[3]);
    }
    else if (weights.ordering_mode == 2)
    {
        // mirror y
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", prim_pts_idx[1], new_edge_pts[0], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", new_edge_pts[1], new_face_pt, new_edge_pts[2], prim_pts_idx[2]);
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[2], prim_pts_idx[3]);
    }
    else if (weights.ordering_mode == 3)
    {
        // mirror x
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", new_edge_pts[0], prim_pts_idx[1], new_edge_pts[1], new_face_pt);
        addprim(0, "poly", new_edge_pts[2], prim_pts_idx[2], new_edge_pts[1], new_face_pt);
        addprim(0, "poly", prim_pts_idx[3], new_edge_pts[2], new_face_pt, new_edge_pts[3]);
    }
    else if (weights.ordering_mode == 4)
    {
        // cyclic
        addprim(0, "poly", prim_pts_idx[0], new_edge_pts[0], new_face_pt, new_edge_pts[3]);
        addprim(0, "poly", prim_pts_idx[1], new_edge_pts[1], new_face_pt, new_edge_pts[0]);
        addprim(0, "poly", prim_pts_idx[2], new_edge_pts[2], new_face_pt, new_edge_pts[1]);
        addprim(0, "poly", prim_pts_idx[3], new_edge_pts[3], new_face_pt, new_edge_pts[2]);
    }
    else if (weights.ordering_mode == 5)
    {
        // opposite
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[0], prim_pts_idx[0]);
        addprim(0, "poly", new_edge_pts[1], new_face_pt, new_edge_pts[0], prim_pts_idx[1]);
        addprim(0, "poly", new_edge_pts[1], new_face_pt, new_edge_pts[2], prim_pts_idx[2]);
        addprim(0, "poly", new_edge_pts[3], new_face_pt, new_edge_pts[2], prim_pts_idx[3]);
    }
    */
}

#endif
