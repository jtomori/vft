import hou
import toolutils
import os
import logging

# logging conffig
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

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

# returns a node which belongs to "vft generator" list
def getOutputGenerator(node):
    generator_nodes = set( ["vft_generate_points"] )
    all_children_nodes = outputChildren(node)
    out = None

    for node in all_children_nodes:
        if node.type().name() in generator_nodes:
            out = node
            break
    
    return out

# find all descending nodes
def outputChildren(node):
    children = list( node.outputs() )
    for node in children:
        new_children = node.outputs()
        if len(new_children) == 0:
            break
        else:
            for child in new_children:
                children.append( child )
                outputChildren(child)
    
    return children

# returns a list of CL statements with fractal function calls from a list of fractal nodes
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
        stack.append(name + default_args)
    
    return stack

# helper func
def clStatementsToString(statements):
    return ";\n".join(statements) + ";\n"

# class that will generate fractal generation CL code that Houdini will read from a string parameter and will execute
class GenerateKernel(object):
    def __init__(self):
        self.vft_root_path = self.getVftRootFromPath( hou.getenv("HOUDINI_PATH") )
        self.vft_kernels_path = os.path.join(self.vft_root_path, "ocl/vft_kernels.cl")
        self.vft_kernels = None # can be initialized by loadKernelsFile member func
        self.vft_kernels_parsed = None # will be filled in by parseKernelsFile member func
    
    # this might not work on Windows
    # extracts path to VFT from os-style paths string
    def getVftRootFromPath(self, path):
        paths = path.split(":")
        
        # this will need to be changed if git repository name changes
        pattern = os.sep + "raymarching" + os.sep + "houdini"

        # find pattern in list of paths
        vft_root = ""
        for path in paths:
            if pattern in path:
                vft_root = path
                break
        
        return vft_root
    
    # loads vft_kernels.cl file into member variable
    def loadKernelsFile(self):
        log.debug("Kernels file loaded")
        with open(self.vft_kernels_path, 'r') as file:
            self.vft_kernels = file.read()

        self.vft_kernels_parsed = self.vft_kernels
        
    
    # parses vft_kernels.cl file and replaces PY_* macros and saves it into member varible
    def parseKernelsFile(self, root_node):
        log.debug("Kernels file parsed")
        # generate fractal stack
        fractals_stack_token = "#define PY_FRACTAL_STACK"

        fractals_stack_cl_code = generateClFractalStack( getInputFractalNodes(root_node) )
        fractals_stack_cl_code = clStatementsToString(fractals_stack_cl_code)
        fractals_stack_cl_code = fractals_stack_token + "\n\n" + fractals_stack_cl_code

        self.vft_kernels_parsed = self.vft_kernels_parsed.replace(fractals_stack_token, fractals_stack_cl_code)

# this func should be used in OpenCL kernelcode parameter Python expression
def fillKernelCodeParm():
    log.debug("OpenCL kernel parm evaluated")
    me = hou.pwd()
    kernel = GenerateKernel()
    kernel.loadKernelsFile()
    kernel.parseKernelsFile( me.parent() )

    return kernel.vft_kernels_parsed

# this func should be called by HDA on On Input Changed event
def nodeRecook(**kwargs):
    log.debug("Node recooked")
    #node.cook(force=True)