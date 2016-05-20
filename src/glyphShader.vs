#version 330

in vec2 glyphVertex;
in vec2 glyphMeshVertex;

in vec2 dataVector;
uniform mat4 transform;
out vec4 vColor;

uniform float glyphScale;

void main() {
    //gl_Position = transform * vec4(PartVertex, 0.0, 1);
    vec2 meshCoord;
    float theta;
    if (abs(dataVector.x)<1E-8) { // avoid dividing by 0 (dividing by -0.0 was causing wrong sense)
        theta = acos(-1.0)/2.0 * dataVector.y/abs(dataVector.y);
    } else {
        theta = atan(dataVector.y/dataVector.x);
        if (dataVector.x<0.0) {
            theta = theta+acos(-1.0); // i.e. that+pi
        }
    }
    float norm = sqrt(dataVector.x*dataVector.x + dataVector.y*dataVector.y);
    meshCoord.x = ( glyphMeshVertex.x * cos(theta) - glyphMeshVertex.y*sin(theta) ) * glyphScale * norm;
    meshCoord.y = ( glyphMeshVertex.x * sin(theta) + glyphMeshVertex.y*cos(theta) ) * glyphScale * norm;
    
    
    vec4 Pos = ( vec4(glyphVertex, 0.0, 1.0) + vec4( meshCoord, 0.0,  0.0) );
    
    
    
    
    gl_Position = transform*Pos;
    
        
    vColor = vec4(1.0,0.9,0.1,1.0);
    

}