const glsl = x => x;


// init
var canvas = document.getElementById("main");
var gl = canvas.getContext('webgl');

// canvas size
// Lookup the size the browser is displaying the canvas.
var displayWidth  = canvas.clientWidth;
var displayHeight = canvas.clientHeight;

// Check if the canvas is not the same size.
if (canvas.width  != displayWidth || canvas.height != displayHeight) {
    // Make the canvas the same size
    canvas.width  = displayWidth;
    canvas.height = displayHeight;
}

/**
 * Shaders
 */

// Utility to fail loudly on shader compilation failure
function compileShader(shaderSource, shaderType) {
    var shader = gl.createShader(shaderType);
    gl.shaderSource(shader, shaderSource);
    gl.compileShader(shader);

    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        throw "Shader compile failed with: " + gl.getShaderInfoLog(shader);
    }

    return shader;
}

var vertexShader = compileShader(glsl`
    attribute vec2 position;

    void main() {
        // position specifies only x and y.
        // We set z to be 0.0, and w to be 1.0
        gl_Position = vec4(position, 0.0, 1.0);
    }
`, gl.VERTEX_SHADER);

var fragmentShader = compileShader(glsl`
    precision highp float;
    uniform vec2 res;

    void main(){
        // Draw every pixel red
        gl_FragColor = vec4(gl_FragCoord.x/res.x, gl_FragCoord.y/res.y, 0.0, 1.0);
    }
`, gl.FRAGMENT_SHADER);

var program = gl.createProgram();
gl.attachShader(program, vertexShader);
gl.attachShader(program, fragmentShader);
gl.linkProgram(program);
gl.useProgram(program);

// setup uniforms
var resLocation = gl.getUniformLocation(program, "res");
gl.uniform2f(resLocation, displayWidth, displayHeight);

/**
 * Geometry setup
 */

// Set up 4 vertices, which we'll draw as a rectangle
// via 2 triangles
//
//   A---C
//   |  /|
//   | / |
//   |/  |
//   B---D
//
// We order them like so, so that when we draw with
// gl.TRIANGLE_STRIP, we draw triangle ABC and BCD.
var vertexData = new Float32Array([
    -1.0,  1.0, // top left
    -1.0, -1.0, // bottom left
     1.0,  1.0, // top right
     1.0, -1.0, // bottom right
]);
var vertexDataBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, vertexDataBuffer);
gl.bufferData(gl.ARRAY_BUFFER, vertexData, gl.STATIC_DRAW);

/**
 * Attribute setup
 */

// Utility to complain loudly if we fail to find the attribute

function getAttribLocation(program, name) {
    var attributeLocation = gl.getAttribLocation(program, name);
    if (attributeLocation === -1) {
        throw 'Can not find attribute ' + name + '.';
    }
    return attributeLocation;
}

// To make the geometry information available in the shader as attributes, we
// need to tell WebGL what the layout of our data in the vertex buffer is.
var positionHandle = getAttribLocation(program, 'position');
gl.enableVertexAttribArray(positionHandle);
gl.vertexAttribPointer(positionHandle,
                       2, // position is a vec2
                       gl.FLOAT, // each component is a float
                       gl.FALSE, // don't normalize values
                       2 * 4, // two 4 byte float components per vertex
                       0 // offset into each span of vertex data
                       );

/**
 * Draw
 */

// Render the 4 vertices specified above (starting at index 0)
// in TRIANGLE_STRIP mode.
gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
