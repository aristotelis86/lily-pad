/**********************************************************************
      Spring class: Creates the springs of a flexible 
      structure to simulate its internal dynamics. 

Example code:
int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution
int N = 2;
PVector gravity = new PVector(0,0);

Window view; 
ControlPoint [] cpoints;
Spring spring;
CollisionHandler collider;

void settings(){
    size(600, 600);
}

void setup() {
  
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  
  cpoints = new ControlPoint[N];
  cpoints[0] = new ControlPoint( new PVector(nx/3.,ny/2.), 5,  10, view );
  cpoints[1] = new ControlPoint( new PVector(2*nx/3.,ny/2.), 5,  10, view );
  
  spring = new Spring(cpoints[0], cpoints[1], nx/8., 5, 0.5, 1, view );
  
  collider = new CollisionHandler( cpoints );
} // end of setup

void draw() {
  background(185);
  for (ControlPoint cp : cpoints) cp.clearForce();
  spring.ApplyAllForces();
  for (ControlPoint cp : cpoints) {
    cp.ApplyForce( gravity );
    cp.updateAlt( 0.1 );
    cp.updateAlt2( 0.1 );
  }
  collider.HandleCollisions();
  for (ControlPoint cp : cpoints) cp.display();
  spring.display();
}
**********************************************************************/

class Spring {
  //========== Attributes - Physical ============//
  float stiffness; // stiffness of spring
  float restLength; // resting length 
  float damping; // damping for simulating the presence of dashpot
  ControlPoint p1, p2; // particle that it is connected to
  
  // For display purposes
  Window myWindow;
  float thick; // default thickness = 1
  color c; // default random coloring
  
  //=============== Constructor =================//
  Spring( ControlPoint a, ControlPoint b, float r, float s, float d, float th_, Window w ) {
    p1 = a;
    p2 = b;
    
    stiffness = s;
    restLength = r;
    damping = d;
    
    myWindow = w;
    thick = th_;
    c = #FF9900; 
  }
  Spring( ControlPoint a, ControlPoint b, float r, float s, float d, Window w) { 
    this( a, b, r, s, d, 1, w);
  }
  Spring( ControlPoint a, ControlPoint b, Window w) { 
    this( a, b, 1, 1, 1, 1, w);
  }
  
  //=================== Methods ================//
  
  // Display
  void display(){
    strokeWeight(thick*myWindow.x.r);
    stroke(c);
    line(myWindow.px(p1.position.x), myWindow.py(p1.position.y), myWindow.px(p2.position.x), myWindow.py(p2.position.y));
  }
  
  // Set thickness
  void matchThickness() { thick = 0.5*(p1.diameter+p2.diameter); }
  
  // Apply Forces on connected particles
  void ApplyAllForces() {
    // apply force due to spring
    PVector springDir = PVector.sub(p1.position, p2.position);
    
    float stretch = springDir.mag();
    stretch -= restLength;
    
    springDir.normalize();
    PVector Tension = PVector.mult(springDir,-stiffness * stretch);
    p1.force.add(Tension);
    Tension.mult(-1);
    p2.force.add(Tension);
    
    // apply force due to dashpot
    PVector RelatVel = PVector.sub(p1.velocity, p2.velocity);
    
    float DampMag = PVector.dot(RelatVel, springDir);
    DampMag = (-1)*DampMag*damping;
    
    PVector DampVec = PVector.mult(springDir, DampMag);
    
    p1.force.add(DampVec);
    DampVec.mult(-1);
    p2.force.add(DampVec);
  }
  
  // Get the stretch of the spring
  float getStretch() {
    PVector SpringDir = PVector.sub(p1.position, p2.position);
    float stretch = SpringDir.mag();
    stretch -= restLength;
    
    return stretch;
  }
  
  // Assign different stiffness from the one constructed
  void updateStiffness(float kk) {
    stiffness = kk;
  }
  
  // Update the damping coefficient if needed
  void UpdateDamping(float dd) {
    damping = dd;
  }
  
  // Export info on the motion of the control points
  void dampInfo( PrintWriter outFile ) {
    outFile.println(getStretch());
  }
  
} // end of Spring class