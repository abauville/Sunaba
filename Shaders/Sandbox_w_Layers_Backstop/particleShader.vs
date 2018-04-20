#version 330

in vec2 PartVertex;
in vec3 PartMeshVertex;
in float PartData;
in float PartPassiveData;
uniform mat4 transform;
out vec4 vColor;
uniform float size;

uniform float type;
uniform vec2 colorScale;


void main() {
    //gl_Position = transform * vec4(PartVertex, 0.0, 1);
    vec4 Pos = ( vec4(PartVertex, 0.0, 1) + vec4(PartMeshVertex, 0) );
    gl_Position = transform*Pos;
    

    if (type==0) {
        // Phase + passive visualization
        if (PartData == 0) {
            vColor = vec4(1.0,1.0,1.0,1.0);
        }
        else if (PartData == 1) {
            //vColor = vec4(0.0,1.0,0.5,1.0);
            vColor = vec4(0.0,.8,0.5,1.0);
            
            if  (PartPassiveData==0.0) {
                vColor = vec4(0.0,.8,0.5,1.0);
            }
            else if (PartPassiveData==1.0) {
                vColor = vec4(0.0,.8*.7,0.5*.7,1.0);
            }
            else if (PartPassiveData==2.0) {
                vColor = vec4(1.0*.75,.9*.75,0.6*.75,1.0);
            }    
            else if (PartPassiveData==3.0) {
                vColor = vec4(1.0,.9,0.6,1.0);
            }
            else if (PartPassiveData==4.0) {
                vColor = vec4(.85,.9,.75,1.0);
            }   
            //vColor = vec4(0.0,.5,0.25,.5);
            //vColor = vec4(0.65,.65,.65,1.0);
        }
        else if (PartData == 2) {
            vColor = vec4(0.1,0.1,.1,1.0);
        }
        else if (PartData == 3) {
            if  (PartPassiveData==0.0) {
                vColor = vec4(0.5,0.4,.8,1.0);
            }
            else if (PartPassiveData==1.0) {
                vColor = vec4(0.5*.7,0.4*.7,.8*.7,1.0);
            }
            else if (PartPassiveData==2.0) {
                vColor = vec4(0.6*.75,1.0*.75,.9*.75,1.0);
            }    
            else if (PartPassiveData==3.0) {
                vColor = vec4(0.6,1.0,.9,1.0);
            }
            else if (PartPassiveData==4.0) {
                vColor = vec4(.75,.85,.9,1.0);
            }   
        }
        else if (PartData == 4) {
            vColor = vec4(0.0,1.0,1.0,1.0);
        }
        else if (PartData == 5) {
            vColor = vec4(1.0,0.0,1.0,1.0);
        }
        else {
            vColor = vec4(1.0,1.0,1.0,1.0);
        }
        
        
        
   } else {
       
       
        // data visualization
        //float pU;
        float ca1, ca2, cb1, cb2, cc1, cc2, cd1, cd2; // limit of caxis for each color
        float R, G, B;
        
       // U = PartData;
        
        ca1 = 0.0;    ca2 =  1.0*colorScale[1];
        cb1 = ca2;    cb2 =  2*ca2;
        cc1 = cb2;    cc2 = 2*cb2; //  with cap at 3*ca1 so that the last color is only 0.5 max
        
        if (PartData>0) {
            R = (PartData-ca1)/(ca2-ca1);
            G = (PartData-cb1)/(cb2-cb1);
            B = (PartData-cc1)/(cc2-cc1);
            
            if (R>1.0){
                R = 1.0;
            }
            if (G>1.0){
                G = 1.0;
            }
            if (B>0.5){
                B = 0.5;
            }
            
        } else {
            B = (PartData+ca1)/(-(ca2-ca1));
            G = (PartData+cb1)/(-(cb2-cb1));
            R = (PartData+cc1)/(-(cc2-cc1));
            
            if (R>0.5){
                R = 0.5;
            }
            if (G>1.0){
                G = 1.0;
            }
            if (B>1.0){
                B = 1.0;
            }
        }
        
        
        if (R<0.0){
            R = 0.0;
        }
        if (G<0.0){
            G = 0.0;
        }
        if (B<0.0){
            B = 0.0;
        }
        
        
        //R = 1-R;
        //B = 1-B;
        //G = 1-G;
        
        // Send color to the fragment shader
        //color.xyz = vec3(R, G, B);
        // color.xyz = vec3(1.0, 1.0, 1.0);
        
        //set every drawn pixel to white
        vColor = vec4(R,G,B,1.0);
      // vColor = vec3(1.0,1.0,1.0);
       

    }
}
