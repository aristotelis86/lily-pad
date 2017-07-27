//=========================== Spring Class ========================//
//**** Extended from original to include dashpot functionality ****//

//-------------------- Example Code for Testing -------------------//
//Window view; // window in h-units
//int m = 40; // x-dir
//int n = 40; // y-dir

//ControlPoint cpoint1, cpoint2;
//Spring spring;

//void setup() {
//  size(500, 500);    
//  frameRate(10);
//  view = new Window(m, n);
  
//  cpoint1 = new ControlPoint(10., 10., 20, view);
//  cpoint2 = new ControlPoint(10., 20., 100, view);
  
//  spring = new Spring( cpoint1, cpoint2, 15, 50, 15, view);
//  spring.thick = view.pdx(.5);
  
//} // end of setup


//void draw() {
//  background(3);
  
//  spring.display();
//  cpoint1.display();
//  cpoint2.display();
  
//  spring.applyAllForces();
  
//  cpoint1.update(0.1);
//  cpoint2.update(0.1);
  
//  println(cpoint1.distance(cpoint2));
//  //noLoop();
//}


class Spring {
  //========== Attributes - Physical ============//
  float stiffness; // stiffness of spring
  float restLength; // resting length 
  float damping; // damping for simulating the presence of dashpot
  ControlPoint p1, p2; // particle that it is connected to
  
  // For display purposes
  Window myWindow;
  float thick; // default thickness = 1
  color c;
  
  //=============== Constructor =================//
  Spring( ControlPoint a, ControlPoint b, float r, float s, float d, Window w) {
    p1 = a;
    p2 = b;
    
    stiffness = s;
    restLength = r;
    damping = d;
    
    myWindow = w;
    thick = myWindow.pdx(1);
    c = color(random(1,255), random(1,255), random(1,255));
    
  }
  
  
  //=================== Methods ================//
  
  // Display
  void display(){
    strokeWeight(thick);
    stroke(c);
    line(myWindow.px(p1.position.x), myWindow.py(p1.position.y), myWindow.px(p2.position.x), myWindow.py(p2.position.y));
  }
  
  // Apply Forces on connected particles
  void applyAllForces() {
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
  
} // end of Spring class