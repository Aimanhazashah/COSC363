/*========================================================================
* COSC 363  Computer Graphics (2018)
* Ray tracer
* See Lab07.pdf for details.
* 
*=========================================================================
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "SceneObject.h"
#include "Ray.h"
#include "TextureBMP.h"
#include "Plane.h"
#include <GL/glut.h>

using namespace std;

const float WIDTH = 20.0;
const float HEIGHT = 20.0;
const float EDIST = 25.0;
const int NUMDIV = 500;
const int MAX_STEPS = 5;
const float XMIN = -WIDTH * 0.5;
const float XMAX =  WIDTH * 0.5;
const float YMIN = -HEIGHT * 0.5;
const float YMAX =  HEIGHT * 0.5;

vector<SceneObject*> sceneObjects;  //A global list containing pointers to objects in the scene


//---The most important function in a ray tracer! ----------------------------------
//   Computes the colour value obtained by tracing a ray and finding its
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{

    glm::vec3 backgroundCol(0);
    glm::vec3 light(10, 40, -3);
    glm::vec3 ambientCol(0.2);   //Ambient color of light

    glm::vec3 colorSum; //Used to store the accumulated colour values

    ray.closestPt(sceneObjects);        //Compute the closest point of intersetion of objects with the ray

    if(ray.xindex == -1) return backgroundCol;      //If there is no intersection return background colour

    glm::vec3 normalVector = sceneObjects[ray.xindex]->normal(ray.xpt); //Normal vector on the sphere at the point of intersection
    glm::vec3 lightVector = light - ray.xpt; //The vector from the point of intersection towards the light source.

    float lightDist = glm::length(lightVector); //Distance between the light and the sphere

    glm::vec3 normalizedLight= glm::normalize(lightVector); //Normalize the light vector

    //Creating a shadow ray from the point of intersection with the ray and object in the scene.
    Ray shadow(ray.xpt, normalizedLight);
    shadow.closestPt(sceneObjects);

    glm::vec3 reflVector = glm::reflect(-lightVector, normalVector); //Reflection vector

    glm::vec3 normalizedReflVector = glm::normalize(reflVector);

    glm::vec3 viewVector = ray.pt- ray.xpt; // view vector = origin - point of intersection

    glm::vec3 normalizedViewVector = glm::normalize(viewVector); // view vector is

    float Dotn2 = glm::dot(normalizedReflVector, normalizedViewVector);

    glm::vec3 specular = glm::vec3(1.f, 1.f, 1.f)*pow((max(Dotn2, 0.f)), 4);

    glm::vec3 normalizeNormal = glm::normalize(normalVector);
    float dotn1 = glm::dot(normalizedLight, normalizeNormal); //This gives the dot product of l.n
    glm::vec3 materialCol = sceneObjects[ray.xindex]->getColor(); //else return object's colour

    if (dotn1 <= 0 || (shadow.xindex > -1 && shadow.xdist < lightDist)){
        colorSum = ambientCol * materialCol;
    } else {
        colorSum = ambientCol * materialCol + (max(dotn1, 0.f) * materialCol) + specular;
    }

    if(ray.xindex == 0 && step < MAX_STEPS)
    {
        glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
        Ray reflectedRay(ray.xpt, reflectedDir);
        glm::vec3 reflectedCol = trace(reflectedRay, step+1); //Recursion!
        colorSum = colorSum + (0.8f*reflectedCol);
    }

    return colorSum;
}

//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
    float xp, yp;  //grid point
    float cellX = (XMAX-XMIN)/NUMDIV;  //cell width
    float cellY = (YMAX-YMIN)/NUMDIV;  //cell height

    glm::vec3 eye(0., 0., 0.);  //The eye position (source of primary rays) is the origin

    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glBegin(GL_QUADS);  //Each cell is a quad.

    for(int i = 0; i < NUMDIV; i++)     //For each grid point xp, yp
    {
        xp = XMIN + i*cellX;
        for(int j = 0; j < NUMDIV; j++)
        {
            yp = YMIN + j*cellY;

            glm::vec3 dir(xp+0.5*cellX, yp+0.5*cellY, -EDIST);  //direction of the primary ray

            Ray ray = Ray(eye, dir);        //Create a ray originating from the camera in the direction 'dir'
            ray.normalize();                //Normalize the direction of the ray to a unit vector
            glm::vec3 col = trace (ray, 1); //Trace the primary ray and get the colour value

            glColor3f(col.r, col.g, col.b);
            glVertex2f(xp, yp);             //Draw each cell with its color value
            glVertex2f(xp+cellX, yp);
            glVertex2f(xp+cellX, yp+cellY);
            glVertex2f(xp, yp+cellY);
        }
    }

    glEnd();
    glFlush();
}

void drawCube(float x, float y, float z, float l, float w, float h, glm::vec3 color) //Fix this shit 
{   
    glm::vec3 A = glm::vec3(x,y,z);
    glm::vec3 B = glm::vec3(x+l,y,z);
    glm::vec3 C = glm::vec3(x+l,y+h,z);
    glm::vec3 D = glm::vec3(x,y+h,z);
    glm::vec3 E = glm::vec3(x+l,y,z-w);
    glm::vec3 F = glm::vec3(x+l,y+h,z-w);
    glm::vec3 G = glm::vec3(x,y+h,z-w);
    glm::vec3 H = glm::vec3(x,y,z-w);
    
    Plane *plane1 = new Plane(A,B,C,D,color);
                              
    Plane *plane2 = new Plane(B,E,F,C,color);
                              
    Plane *plane3 = new Plane(E,H,G,F,color);
    
    Plane *plane4 = new Plane(D,G,H,A,color);
    
    Plane *plane5 = new Plane(D,C,F,G,color);
    
    Plane *plane6 = new Plane(H,E,B,A,color);
                                                      
    sceneObjects.push_back(plane1);
    sceneObjects.push_back(plane2);
    sceneObjects.push_back(plane3);
    sceneObjects.push_back(plane4);
    sceneObjects.push_back(plane5);
    sceneObjects.push_back(plane6);
    
}




//---This function initializes the scene -------------------------------------------
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize()
{
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    //texture = TextureBMP("Earth.bmp") //texture
    glClearColor(0, 0, 0, 1);

    Plane *plane = new Plane (glm::vec3(-30., -20, -30), //Point A
                                glm::vec3(30., -20, -30), //Point B
                                glm::vec3(30., -20, -100), //Point C
                                glm::vec3(-30., -20, -100), //Point D
                                glm::vec3(0.5, 1, 1)); //Colour
                                
    //-- Create a pointer to a sphere object
    Sphere *sphere1 = new Sphere(glm::vec3(-5.0, -5.0, -90.0), 15.0, glm::vec3(0, 0, 1));
    Sphere *sphere2 = new Sphere(glm::vec3(-8.5, -3.5, -65.0), 2.0, glm::vec3(1, 0, 0));
    Sphere *sphere3 = new Sphere(glm::vec3(1.5, -10.0, -65.0), 3.0, glm::vec3(0, 1, 0));
    //--Add the above to the list of scene objects.
    sceneObjects.push_back(sphere1);
    sceneObjects.push_back(sphere2);
    sceneObjects.push_back(sphere3);
    sceneObjects.push_back(plane);
    drawCube(8, -10.0, -60.0,6,10,6,glm::vec3(1, 0.7529, 0.7961));//handling the cube


}



int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracer");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
