<Please submit this file with your solution.>

2020 Spring CSCI 520, Assignment 1

<Richard Chien>

================
<Coding Environment>

Visual Studio 2017  | x86 | Debug

<Description of what you have accomplished>

1. Mass spring system.(structural,shear,bend spring)---I use Hook's law for creating the spring force and dampling force to slow down the motion. Also, I implement these force on the network of spring(structural,shear,bend spring). All are followed from the slides from helper slides.

2. Incorporate a bounding box, including collision detection and response--- I set up a CollisionDetection function to detect the collision. Once collide, the penalty method and collision spring are applied. 

3. Implement an arbitrary non-homogeneous time-independent force field, with appropriate force-field interpolation(ex. external.w)--- set the IsOtherForce function to detect whether has force applied. 

If it is, I calculate force at each point by interpolating the eight force nodes around the cube. Apply that force in each point in cube.

4. Able to produce JPEG frames and do the animation.

Related function:
bool CollisionDetection(struct point pos,double boxsize,int *result);
bool IsOtherForce(struct world * jello);
Vector PenaltyForce(struct world * jello,int i,int j, int k);
Vector ExternalForce(struct world * jello,int i,int j, int k);
Vector ForceInSpring(struct world * jello,double RLength,struct point pos1, struct point pos2, struct point V1,struct point V2);

<Also, explain any extra credit that you have implemented.>

1. If the world file has external force, you can click 'f' to deactivate it or active it again.

2. Implement collision detection with inclined plane
The related world file is incline.w(has an inclined plane in X+Y=0) or inclinez.w(z=-1.5)
First check every point whether they are on the same side of inclined plane. If the points in cube cross the plain then give it a force to push it back. And that force is produced by calculating the L vector and same method as incorporating a bounding box. 

3. Drag the mouse to give the force to all the vertices in cube.
Press the left and drag the mouse in any direction and release will give the cube an external force. I also printf the mouse location from start to end. To do this, by calculating the length from start point to end point. And give a coefficient to modify the force, then apply that force to every points in cube.


4. Set up texture mapping to the cube. 
If press "t" and the cube is at triangle rendering mode. It will apply the dice texture to the cube. 
I accomplish this by load six dice images first, then use glGenTextures(1, &eachTex); glBindTexture(GL_TEXTURE_2D, eachTex),etc. init function to load texture, then apply into the cube when press the ¨t¨. 

f: activate or deactivate the external force
t: display checkerboard texture or hide it

Related function:
bool IsExternalForce(struct world * jello);
Vector InclinedForce(struct world * jello,int i,int j, int k);
Vector MouseForce(struct world * jello);
void showTexture(struct world * jello);

<Animation>
In 300 pictures, first I give the cube velocity to collide the invisable inclined plain(X+Y=0). Then, when the cube is pushed back, I drag the cube from up to down.
Next, I press "t" in triangle rendering to show my texture. 

