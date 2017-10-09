/**********************************************************************
      ControlPoint class: Creates the control points of 
      a flexible structure and handles their behaviour
      during the runs. 

Example code:
int nx = (int)pow(2,4); // x-dir
int ny = (int)pow(2,4); // y-dir

int N = 5; // number of control points
PVector gravity = new PVector(2,4);

Window view; 
ControlPoint [] cpoints = new ControlPoint[N];
CollisionHandler collide;

void settings(){
    size(600, 600);
}

void setup() {
  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  for (int i=0; i<N; i++) {
    cpoints[i] = new ControlPoint( new PVector(random(view.x.inE),random(view.y.inE)), 5,  10, view );
  }
  collide = new CollisionHandler( cpoints );
}

void draw() {
  background(185);
  for (ControlPoint cp : cpoints) {
    cp.clearForce();
    cp.ApplyForce( gravity );
    cp.updateAlt( 0.05 );
    cp.updateAlt2( 0.05 );
  }
   
  collide.HandleCollisions();
  for (ControlPoint cp : cpoints) cp.display();
}
**********************************************************************/

class ControlPoint {
  //================= Attributes ====================//
  
  PVector position; // current position
  PVector positionOld; // position before updating
  PVector velocity; // current velocity
  PVector velocityOld; // velocity before updating
  PVector acceleration; // current acceleration
  PVector accelerationOld; // acceleration to begin the update
  PVector force; // force acting on the point-mass
  PVector impForce; // impact force tracking
  float mass; // mass of the point
  boolean xfixed; // fix the particle at its y-axis
  boolean yfixed; // fix the particle at its x-ayis
  
  Window myWindow; // viewing window
  color c; // for displaying
  float diameter; // for collisions model (mostly)
  
  //================= Constructor ====================//
  ControlPoint(PVector position_, float m,  float thk, Window myWindow_) {
    position = position_;
    velocity = new PVector(0, 0);
    acceleration = new PVector(0, 0);
    force = new PVector(0, 0);
    mass = m;
    
    xfixed = false;
    yfixed = false;
    
    myWindow = myWindow_;
    c = color(random(1,255), random(1,255), random(1,255));
    diameter = thk;
    
    positionOld = position.copy();
    velocityOld = velocity.copy();
    accelerationOld = acceleration.copy();
  }
  
  ControlPoint(PVector position_, float m, Window myWindow_) {this(position_, m,  1, myWindow_);}
  
  
  //================= Methods =====================//
  
  // Display (for testing)
  void display() {
    noStroke();
    fill(c);
    ellipse(myWindow.px(position.x), myWindow.py(position.y), diameter*myWindow.x.r, diameter*myWindow.y.r);
  }
  
  // Display impact (for testing)
  void impDisplay() {
    noStroke();
    fill(255, 0, 0);
    ellipse(myWindow.px(position.x), myWindow.py(position.y), 0.8*diameter*myWindow.x.r, 0.8*diameter*myWindow.y.r);
  }
  
  // Update the position
  void UpdatePosition(float x, float y) {
    position.x = x;
    position.y = y;
  }
  
  // Update the velocity
  void UpdateVelocity(float x, float y) {
    velocity.x = x;
    velocity.y = y;
  }
  
  // Store old state
  void StoreOld() {
    positionOld = position.copy();
    velocityOld = velocity.copy();
    accelerationOld = acceleration.copy();
  }
  
  // For testing...
  void randomWalk() {
    this.StoreOld();
    float ex, ey, vx, vy, x, y;
    ex = random(1);
    ey = random(1);
    
    if (ex<.05) vx = -1*velocity.x;
    else vx = velocity.x;
    if (ey<.05) vy = -1*velocity.y;
    else vy = velocity.y;
    
    this.UpdateVelocity(vx,vy);
    x = positionOld.x + velocity.x;
    y = positionOld.y + velocity.y;
    this.UpdatePosition(x,y);
  }
  
  // Clear any forces acting on the particle
  void clearForce() { force.mult(0); }
  
  // Accumulate all the forces acting on the particle
  void ApplyForce(PVector FF) { force.add( FF ); }
  
  // Find the acceleration due to forces
  void calculateAcceleration() {
    PVector accel = force.copy();
    acceleration = accel.div(mass);
  }
  
  // Make the particle free of constraints
  void makeFree() {
    xfixed = false;
    yfixed = false;
  }
  void makeFreex() { xfixed = false; }
  void makeFreey() { yfixed = false; }
  
  // Constrain the particle at its location
  void makeFixed() {
    xfixed = true;
    yfixed = true;
  }
  void makeFixedx() { xfixed = true; }
  void makeFixedy() { yfixed = true; }
  
  boolean isFixedx() { return xfixed; }
  boolean isFixedy() { return yfixed; }
  
  // Update methods based on Predictor-Corrector scheme 
  void update( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      StoreOld();
      float x, y, vx, vy;
      
      if (!this.isFixedx()) {
        x = position.x + t*velocity.x;
        vx = velocity.x + t*acceleration.x;
      }
      else {
        x = position.x;
        vx = velocity.x;
      }
      if (!this.isFixedy()) {
        y = position.y + t*velocity.y;
        vy = velocity.y + t*acceleration.y;
      }
      else {
        y = position.y;
        vy = velocity.y;
      }
      UpdatePosition( x, y );
      UpdateVelocity( vx, vy );
    }
  }
  
  void update2( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      float x, y, vx, vy;
      if (!this.isFixedx()) {
        x = positionOld.x + .5*t*(velocityOld.x + velocity.x);
        vx = velocityOld.x + .5*t*(accelerationOld.x + acceleration.x);
      }
      else {
        x = position.x;
        vx = velocity.x;
      }
      if (!this.isFixedy()) {
        y = positionOld.y + .5*t*(velocityOld.y + velocity.y);
        vy = velocityOld.y + .5*t*(accelerationOld.y + acceleration.y);
      }
      else {
        y = position.y;
        vy = velocity.y;
      }
      UpdatePosition( x, y );
      UpdateVelocity( vx, vy );
    }
  }
  
  // Alternative update methods based on Predictor-Corrector scheme
  // The position is already 2nd order
  void updateAlt( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      StoreOld();
      float x, y, vx, vy;
      
      if (!this.isFixedx()) {
        x = position.x + t*velocity.x + 0.5*acceleration.x*sq(t);
        vx = velocity.x + t*acceleration.x;
      }
      else {
        x = position.x;
        vx = velocity.x;
      }
      if (!this.isFixedy()) {
        y = position.y + t*velocity.y + 0.5*acceleration.y*sq(t);
        vy = velocity.y + t*acceleration.y;
      }
      else {
        y = position.y;
        vy = velocity.y;
      }
      UpdatePosition( x, y );
      UpdateVelocity( vx, vy );
    }
  }
  // only update velocity here...
  void updateAlt2( float t ) {
    if ((!this.isFixedx()) || (!this.isFixedy())) {
      calculateAcceleration();
      float vx, vy;
      if (!this.isFixedx()) vx = velocityOld.x + .5*t*(accelerationOld.x + acceleration.x);
      else vx = velocity.x;
      if (!this.isFixedy()) vy = velocityOld.y + .5*t*(accelerationOld.y + acceleration.y);
      else vy = velocity.y;
      UpdateVelocity( vx, vy );
    }
  }
  
  // Export info on the motion of the control points
  void dampInfo( PrintWriter outFile ) {
    outFile.println(position.x+","+position.y+","+velocity.x+","+velocity.y+","+force.x+","+force.y);
  }
  void dampInfo( PrintWriter outFile, float t ) {
    outFile.println(t+","+position.x+","+position.y+","+velocity.x+","+velocity.y+","+force.x+","+force.y);
  }
  
  // Get the distance between control points
  float distance(ControlPoint other) {
    float d = this.position.dist(other.position);
    return d;
  }
  
  // Rewind the update
  void rewindPosition( float dt ) {
    float x, y;
    x = position.x - dt*(position.x - positionOld.x);
    y = position.y - dt*(position.y - positionOld.y);
    UpdatePosition(x,y);
  }
  
  // Apply impact forces to main force variable
  void applyImpForce() { force.add( impForce ); }
  
} // end of ControlPoint class
