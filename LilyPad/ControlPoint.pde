//======================= ControlPoint Class ======================= //

// ------------------- Example Code for Testing ------------------- //
//Window view; // window in h-units
//int m = 40; // x-dir
//int n = 40; // y-dir

//ControlPoint [] cpoint;

//void setup() {
//  size(500, 500);    
//  frameRate(10);
//  view = new Window(m, n);
  
//  cpoint = new ControlPoint[3];
//  for (int i=0; i<3; i++) cpoint[i] = new ControlPoint(2., 5+i*10., 20+i*30, view);
    
  
//} // end of setup


//void draw() {
//  background(3);
  
//  for (ControlPoint c : cpoint) c.display();
  
//  cpoint[0].applyForce(new PVector(1,0));
//  cpoint[1].applyForce(new PVector(1,0));
//  cpoint[2].applyForce(new PVector(1,0));
  
//  for (ControlPoint c : cpoint) c.update(1);
  
//  //noLoop();
//}



class ControlPoint {
  //================= Attributes ====================//
  
  PVector position; // current position
  PVector positionOld;
  PVector velocity; // current velocity
  PVector velocityOld;
  PVector force; // force acting on the point-mass
  float mass; // mass of the point
  boolean fixed; // fix the particle at its location
  boolean xfixed; // fix the particle at its y-axis
  boolean yfixed; // fix the particle at its x-ayis
  
  Window myWindow; // viewing window
  color c; // for displaying
  float thick;
  
  //================= Constructor ====================//
  ControlPoint(PVector position_, float m,  float thk, Window myWindow_) {
    position = position_;
    velocity = new PVector(0, 0);
    force = new PVector(0, 0);
    mass = m;
    
    fixed = false;
    xfixed = false;
    yfixed = false;
    
    myWindow = myWindow_;
    c = color(random(1,255), random(1,255), random(1,255));
    thick = thk;
    
    positionOld = position.copy();
    velocityOld = velocity.copy();
  }
  
  ControlPoint(PVector position_, float m, Window myWindow_) {this(position_, m,  1, myWindow_);}
  
  
  //================= Methods =====================//
  
  // Display it
  void display() {
    noStroke();
    fill(c);
    ellipse(myWindow.px(position.x), myWindow.py(position.y), myWindow.x.r, myWindow.y.r);
  }
  
  // Clear any forces acting on the particle
  void clearForce() {
    force.mult(0);
  }
  
  // Accumulate all the forces acting on the particle
  void applyForce(PVector FF) {
    force.add(FF);
  }
  
  // Make the particle free of constraints
  void makeFree() {
    fixed = false;
    xfixed = false;
    yfixed = false;
  }
  
  void makeFreex() { xfixed = false; }
  void makeFreey() { yfixed = false; }
  
  // Constrain the particle at its location
  void makeFixed() {
    fixed = true;
    xfixed = true;
    yfixed = true;
  }
  
  // Constrain the particle at its y-axis
  void makeFixedx() { xfixed = true; }
  
  // Constrain the particle at its x-axis
  void makeFixedy() { yfixed = true; }
  
  
  // Get the distance between control points
  float distance(ControlPoint other) {
    float d = this.position.dist(other.position);
    return d;
  }
  
  // For testing...
  void move() {
    PVector randVel = new PVector(random(-1,1), random(-1,1));
    //position.add(randVel);
    position.add(new PVector(1,0));
  }
  
  void update(float t) {
    PVector accel = PVector.div(force, mass);
    
    velocity.add(PVector.mult(accel,t));
    position.add(PVector.mult(velocity,t));
    
    this.clearForce();
  }
  
  void UpdatePosition(float x, float y) {
    position.x = x;
    position.y = y;
  }
  
  void UpdateVelocity(float x, float y) {
    velocity.x = x;
    velocity.y = y;
  }
} // end of ControlPoint class