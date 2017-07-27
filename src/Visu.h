/*
 * Visu.h
 *
 *  Created on: Jul 27, 2017
 *      Author: abauville
 */

#ifndef VISU_H_
#define VISU_H_

#include "stokes.h"

// Visualization
// ========================
typedef enum { Blank,
               Viscosity,
               StrainRate,
               Velocity,
               Pressure,
               Density,
               Temperature,
               Stress,
               FluidPressure,
               Permeability,
               Porosity,
               CompactionPressure,
               Phase,
               VxRes,
               VyRes,
               PRes,
               PfRes,
               PcRes,
               TRes,
               VelocityDiv,
               SIIOvYield,
               PeOvYield,
               Khi,
               Khib,
               Strain,
               Vorticity,
               POvPlitho,
               EffectiveViscosity,
               ShearModulus } VisuType;
typedef enum { PartPhase,
               PartTemp,
               PartSigma_xx,
               PartSigma_xy,
               PartDeltaP,
               PartPorosity } ParticleVisuType;
typedef enum { StokesVelocity,
               DarcyGradient,
               DeviatoricStressTensor } GlyphType;
typedef enum { Triangle,
               ThinArrow,
               ThickArrow,
               TensorCross } GlyphMeshType;
typedef enum { Nearest,
               Linear } FilterType;

//typedef enum {Visu_Alpha_Phase, Visu_Alpha_Threshold, Visu_Alpha_AbsThreshold} VisuAlphaType;

struct ColorMap
{
    //number     = number
    //type       = colormapType # "automatic would go from min to max values"
    compute colorMapRes;
    compute colorMap;
    compute scale;
    compute center; // centered value (scaled)
    compute max;    // maximum value (scaled) from the center
    bool log10on;
    compute alphaAbsThreshold; // absolute value of the threshold for transparecny (not affected by log10on)
};

struct Visu
{

    GLFWwindow *window;

    int ntri, ntrivert;
    GLuint *elements;
    GLfloat *U;
    GLfloat *vertices;
    GLfloat scale, valueScale, valueShift;
    GLfloat colorScale[2], partColorScale[2];
    GLfloat shift[3], shiftFac[3];
    GLint log10_on;
    GLuint VAO, VBO, EBO;
    GLuint TEX;
    GLuint VAO_part, VBO_part, VBO_partMesh;
    GLuint VAO_glyph, VBO_glyph, VBO_glyphMesh;
    GLuint ShaderProgram, ParticleShaderProgram, ParticleBackgroundShaderProgram, GlyphShaderProgram;

    bool closeAtTheEndOfSimulation;

    /*
	const char* VertexShaderFile;
	const char* FragmentShaderFile;

	const char* ParticleVertexShaderFile;
	const char* ParticleFragmentShaderFile;
	const char* ParticleGeometryShaderFile;
	const char* ParticleBackgroundVertexShaderFile;
	const char* ParticleBackgroundFragmentShaderFile;
	const char* GlyphVertexShaderFile;
	const char* GlyphFragmentShaderFile;
	*/

    char VertexShaderFile[MAX_STRING_LENGTH];
    char FragmentShaderFile[MAX_STRING_LENGTH];

    char ParticleVertexShaderFile[MAX_STRING_LENGTH];
    char ParticleFragmentShaderFile[MAX_STRING_LENGTH];
    char ParticleGeometryShaderFile[MAX_STRING_LENGTH];
    char ParticleBackgroundVertexShaderFile[MAX_STRING_LENGTH];
    char ParticleBackgroundFragmentShaderFile[MAX_STRING_LENGTH];
    char GlyphVertexShaderFile[MAX_STRING_LENGTH];
    char GlyphFragmentShaderFile[MAX_STRING_LENGTH];

    VisuType type;
    ParticleVisuType typeParticles;

    GLfloat *particles;
    int nParticles;
    bool showParticles;
    GLfloat *particleMesh;
    int particleMeshRes;
    compute particleMeshSize;
    int nGlyphs;
    int glyphSamplingRateX;
    int glyphSamplingRateY;

    GLfloat glyphScale;

    GLfloat *glyphs;
    GLfloat *glyphMesh;

    // Input variables
    bool mouse1Pressed;
    bool mouse2Pressed;
    double mouse1BeginDrag[2];
    double mouse1EndDrag[2];
    double mouse2BeginDrag[2];
    double mouse2EndDrag[2];
    double mouse2BeginDragShifted[2];

    bool paused;
    bool nonLinItisOver;
    bool initPassivePart;

    GLFWcursor *handCursor;

    unsigned char *imageBuffer; // stores the pixel data from the window to be stored in an image file
    bool writeImages;

    int retinaScale;

    char outputFolder[MAX_STRING_LENGTH];
    char shaderFolder[MAX_STRING_LENGTH];

    bool transparency;
    bool alphaOnValue;

    bool showGlyphs;
    GlyphType glyphType;
    GlyphMeshType glyphMeshType;
    int nGlyphMeshVert;

    bool update;

    int width, height;

    bool updateGrid;

    FilterType filter;
    compute alphaAbsThreshold;

    ColorMap colorMap[MAX_VISU_TYPE];

};

// Visualization
// =========================
void Visu_Memory_allocate(Visu *Visu, Grid *Grid);
void Visu_Memory_free(Visu *Visu);
void Visu_init(Visu *Visu, Grid *Grid, Particles *Particles, Char *Char, Input *Input);
void Visu_updateVertices(Visu *Visu, Grid *Grid);
void Visu_initWindow(Visu *Visu);
void error_callback(int error, const char *description);
void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);

void Visu_updateInterp_Any_Node2Cell_Local(Visu *Visu, Grid *Grid, compute *CellValue);
void Visu_updateInterp_Any_Node2Cell_Locali(Visu *Visu, Grid *Grid, int *CellValue);
void Visu_StrainRate(Visu *Visu, Grid *Grid, Physics *Physics);
void Visu_updateUniforms(Visu *Visu);
void Visu_velocity(Visu *Visu, Grid *Grid, Physics *Physics);
void VisudivV(Visu *Visu, Grid *Grid, Physics *Physics);
void Visu_stress(Visu *Visu, Grid *Grid, Physics *Physics);
void Visu_SIIOvYield(Visu *Visu, Grid *Grid, Physics *Physics, Numerics *Numerics, MatProps *MatProps);
void Visu_POvPlitho(Visu *Visu, Grid *Grid, Physics *Physics, Numerics *Numerics);
void Visu_PeOvYield(Visu *Visu, Grid *Grid, Physics *Physics, Numerics *Numerics);
void Visu_update(Visu *Visu, Grid *Grid, Physics *Physics, Char *Char, EqSystem *EqStokes, EqSystem *EqThermal, Numbering *NumStokes, Numbering *NumThermal, Numerics *Numerics, MatProps *MatProps);
void Visu_checkInput(Visu *Visu);
void Visu_particles(Visu *Visu, Particles *Particles, Grid *Grid);
void Visu_glyphs(Visu *Visu, Physics *Physics, Grid *Grid, Particles *Particles);
void Visu_particleMesh(Visu *Visu);
void Visu_alphaValue(Visu *Visu, Grid *Grid, Physics *Physics);
void Visu_glyphMesh(Visu *Visu);

void Visu_main(Visu *Visu, Grid *Grid, Physics *Physics, Particles *Particles, Numerics *Numerics, Char *Char, EqSystem *EqStokes, EqSystem *EqThermal, Numbering *NumStokes, Numbering *NumThermal, MatProps *MatProps);
void Visu_residual(Visu *Visu, Grid *Grid, EqSystem *EqSystem, Numbering *Numbering);

#endif /* VISU_H_ */
