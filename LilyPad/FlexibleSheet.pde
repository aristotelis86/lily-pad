//==================== FlexibleSheet Class ==================//



class FlexibleSheet extends LineSegBody {
  //============================= Attributes =================================//
  int numOfpoints;
  int numOfsprings;
  float Mass, segLength, pointMass, stiffness;
  float [] xcurrent, ycurrent, vxcurrent, vycurrent;
  float [] xnew, ynew, vxnew, vynew;
  
  // Dependancy from other objects
  ControlPoint [] cpoints; // stores the control points of the system
  Spring [] springs; // stores the connecting springs
  
  FlexibleSheet( float L_, float M_, float stiffness_, float x, float y, PVector align, Window window ) {
    super(x, y, window); 
    weight = window.pdx(thk);
    
    align.normalize();
    PVector X0;
    for ( int i=0; i < L_; i++ ) {
      X0 = PVector.add(new PVector(x,y), PVector.mult(align, i));
      super.add(X0.x, X0.y);
    }
    super.end();
    
    Mass = M_;
    stiffness = stiffness_;
    numOfpoints = this.coords.size();
    numOfsprings = numOfpoints - 1;
    segLength = 1; // resting length of each spring
    pointMass = Mass/numOfpoints;
    
    cpoints = new ControlPoint[numOfpoints];
    springs = new Spring[numOfsprings];
    
    for (int i = 0; i < numOfpoints; i++) cpoints[i] = new ControlPoint(this.coords.get(i), pointMass, window);
    for (int i = 0; i < numOfsprings; i++) springs[i] = new Spring( cpoints[i], cpoints[i+1], segLength, stiffness, 0, window);
    
    xcurrent = new float[numOfpoints]; xnew = new float[numOfpoints];
    ycurrent = new float[numOfpoints]; ynew = new float[numOfpoints];
    vxcurrent = new float[numOfpoints]; vxnew = new float[numOfpoints];
    vycurrent = new float[numOfpoints]; vynew = new float[numOfpoints];
    
  } // end of Constructor
  
  //============================= Methods =================================//
  boolean unsteady() {return true;}
  
  float velocity( int d, float dt, float x, float y ) {
    int coorSize = super.coords.size();
    float sumWeights = 0;
    float sumVel = 0;
    
    for (int j=0; j<coorSize; j++) {
      PVector r = new PVector(x, y);
      r.sub(coords.get(j));
      float factor = exp((-1)*r.magSq());
      
      if (d==1) {
        sumVel += cpoints[j].velocity.x * factor;
        sumWeights += factor;
      }
      else {
        sumVel += cpoints[j].velocity.y * factor;
        sumWeights += factor;
      }
    }
    if (sumWeights<1e-9) return 0.0;
    else return sumVel/sumWeights;
  }
  
  // compute body reaction to applied force 
  void react (PVector [] extForce, float dt) {
    
  }
  
  
  void move() {
    for (ControlPoint cp : cpoints) {
      cp.velocity = new PVector(0, .2);
      cp.position.add( cp.velocity );
    }
    for ( OrthoNormal o: orth   ) o.translate(0, .2);
    if (n>4) box.translate(0, .2);
  }

} //=========== end of FlexibleSheet class ===============