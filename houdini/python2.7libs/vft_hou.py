import hou
import toolutils
import os
import logging
import time

"""
todo
    * add unique identifier to detail attrib fractal_name
"""

# logging config
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

# returns list of fractal nodes that are connected (upstream) to root node
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

# returns a connected node (downstream) which belongs to "vft generator" list
def getOutputNodeByTypeName(start_node, type_name=""):
    all_children_nodes = outputChildren(start_node)
    out = None

    for node in all_children_nodes:
        if node.type().name() == type_name:
            out = node
            break
    
    return out

# find all descending connected (downstream) nodes
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
    for name in fractal_nodes:
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

        self.vft_kernels = None
        self.vft_kernels_parsed = None
    
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
    def loadKernelsFileToMemberVar(self):
        start_time = time.time()
        with open(self.vft_kernels_path, 'r') as file:
            self.vft_kernels = file.read()


        log.debug("Kernels file loaded from disk in {0:.8f} seconds".format( time.time() - start_time ))
    
    # loads vft_kernels.cl into specified parm object (which should be string) - this function should be called by a button for (re)loading a parm
    def loadKernelsFileToParm(self, parm):
        if self.vft_kernels == None:
            self.loadKernelsFileToMemberVar()

        parm.set(self.vft_kernels)
    
    # loads vft_kernels.cl into member var - either from disk, or parm (if it is loaded there already)
    def loadKernelsFileFromParm(self, parm):
        if parm.eval() == "":
            log.debug("Loading member var from file")
            self.loadKernelsFileToMemberVar()
        else:
            log.debug("Loading member var from node parameter")
            self.vft_kernels = parm.eval()
    
    # parses vft_kernels.cl file and replaces PY_* macros and saves it into member varible
    def parseKernelsFile(self, fractal_nodes):
        start_time = time.time()
        self.vft_kernels_parsed = self.vft_kernels

        # generate fractal stack
        fractals_stack_token = "#define PY_FRACTAL_STACK"

        fractals_stack_cl_code = clStatementsToString( generateClFractalStack(fractal_nodes) )
        fractals_stack_cl_code = fractals_stack_token + "\n\n" + fractals_stack_cl_code

        self.vft_kernels_parsed = self.vft_kernels_parsed.replace(fractals_stack_token, fractals_stack_cl_code)


        log.debug("Kernels file parsed in {0:.8f} seconds".format( time.time() - start_time ))

# this func will do all the parsing and will set up the kernel parm in descendant opencl node
def fillKernelCodePythonSop():
    start_time = time.time()    
    me = hou.pwd()
    geo = me.geometry()
    kernels_parm = me.parm("vft_kernels")

    # find a opencl downstream node
    cl_node = getOutputNodeByTypeName(me, "opencl")

    # init a GenerateKernel object and init member var vft_kernels
    kernel = GenerateKernel()
    kernel.loadKernelsFileFromParm(kernels_parm)

    # get set of incoming fractals
    fractal_nodes = list( geo.findGlobalAttrib("fractal_name").strings() )
    print fractal_nodes

    # do the parsing
    kernel.parseKernelsFile(fractal_nodes)

    # set vft_kernels_parsed to kernelcode parm in an opencl node
    cl_node.parm("kernelcode").set(kernel.vft_kernels_parsed)


    log.debug("Python SOP evaluated in {0:.8f} seconds".format( time.time() - start_time ))    

# this func should be called by HDA on On Input Changed event
def nodeRecook(**kwargs):
    log.debug("Node recooked")
    #node.cook(force=True)