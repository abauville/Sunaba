#version 330 core

in vec3 color;
out vec4 out_Color;

uniform sampler2D tex;

in vec2 TexCoord;

uniform float valueScale;
uniform vec2 colorScale;
uniform int log10_on;
uniform float one_ov_log_of_10;
uniform float valueShift;

void main() {
    
    
    
    
    
    // Final Position of the Vertices
    float pU;
    float U;
    float alpha;
    //    U = 100.0;
    U = texture(tex, TexCoord).r;
    
    // Compute values of color according to the solution
    float ca1, ca2, cb1, cb2, cc1, cc2, cd1, cd2; // limit of caxis for each color
    float R, G, B;
    //ca1 = 0.5;      ca2 = 1;
    //cb1 = 1;        cb2 = 2.5;
    //cc1 = -0.5;     cc2 = -1;
    
    if (log10_on==1) {
        pU = one_ov_log_of_10 * log(U/valueScale);
    }
    else {
        pU = U/valueScale;
    }
    pU += valueShift;
    
    ca1 = 0.0;      ca2 =  1.0*colorScale[1];
    cb1 = ca2;    cb2 =  2*ca2;
    cc1 = cb2;    cc2 = 2*cb2; //  with cap at 3*ca1 so that the last color is only 0.5 max
    
    if (pU>0) {
        R = (pU-ca1)/(ca2-ca1);
        G = (pU-cb1)/(cb2-cb1);
        B = (pU-cc1)/(cc2-cc1);
        
        if (R>1.0){
            R = 1.0;
        }
        if (G>1.0){
            G = 1.0;
        }
        if (B>0.0){
            B = 0.0;
        }
        
    } else {
        B = (pU+ca1)/(-(ca2-ca1));
        G = (pU+cb1)/(-(cb2-cb1));
        R = (pU+cc1)/(-(cc2-cc1));
        
        if (R>0.0){
            R = 0.0;
        }
        if (G>1.0){
            G = 1.0;
        }
        if (B>1.0){
            B = 1.0;
        }
    }
    
    
    
    //if (G<0.0) {
    //    G = (pU+cb1)/(-(cb2-cb1));
    //}
    
    //R = pU;
    
    
    
    
    if (R<0.0){
        R = 0.0;
    }
    if (G<0.0){
        G = 0.0;
    }
    if (B<0.0){
        B = 0.0;
    }
    
    
    
    
    alpha = texture(tex, TexCoord).g;
    
    alpha *= R;
    if (alpha<0) {
        alpha = 0;
    } else if (alpha>1) {
        alpha = 1;
    }
    //alpha = 1.0;
    
    
    
    out_Color = vec4(R,G,B,alpha);
}
