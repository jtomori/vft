import hou
import toolutils

# returns list of fractal nodes that are connected to root node
def getInputFractalNodes(root):
    # node type names of all fractal nodes
    fractals_nodes = set( ["vft_bristorbrotIter", "vft_mandelbulbPower2Iter", "vft_mengerSpongeIter"] )

    all_input_nodes = root.inputAncestors()
    input_fractal_nodes = []

    # find fractal nodes in all input nodes
    for node in all_input_nodes:
        if node.type().name() in fractals_nodes:
            input_fractal_nodes.append(node)

    return input_fractal_nodes

# returns a list of CL fractal function calls
def generateClFractalStack(fractal_nodes):
    default_args = "(Z, de, P_in, log_lin, 1.0f, (float4)(0.0f))"

    fractal_names = []

    # convert houdini node name to fractal function name
    for node in fractal_nodes:
        name = node.type().name()
        name = name.split("_")[-1]
        fractal_names.append(name)
    
    # list which will hold CL fractal funcs calls
    stack = []

    for name in fractal_names:
        stack.append(name + default_args + ";")
    
    return stack