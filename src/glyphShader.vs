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
    
    float theta = atan(dataVector.y/dataVector.x);
    float norm = sqrt(dataVector.x*dataVector.x + dataVector.y*dataVector.y);
    meshCoord.x = ( glyphMeshVertex.x * cos(theta) - glyphMeshVertex.y*sin(theta) ) * glyphScale * norm;
    meshCoord.y = ( glyphMeshVertex.x * sin(theta) + glyphMeshVertex.y*cos(theta) ) * glyphScale * norm;
    
    
    vec4 Pos = ( vec4(glyphVertex, 0.0, 1.0) + vec4( meshCoord, 0.0,  0.0) );
    
    
    
    
    gl_Position = transform*Pos;
    
        
    vColor = vec4(1.0,1.0,1.0,1.0);
    

}