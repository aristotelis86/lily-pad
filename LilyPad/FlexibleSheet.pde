//==================== FlexibleSheet Class ==================//



class FlexibleSheet extends LineSegBody {
  //============================= Attributes =================================//
  int numOfpoints;
  int numOfsprings;
  float Length, Mass, segLength, pointMass, stiffness, damping;
  float [] xcurrent, ycurrent, vxcurrent, vycurrent;
  float [] xnew, ynew, vxnew, vynew;
  float dtmax;
  
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
    
    Length = L_;
    Mass = M_;
    stiffness = stiffness_;
    numOfpoints = this.coords.size();
    numOfsprings = numOfpoints - 1;
    segLength = 1; // resting length of each spring
    pointMass = Mass/numOfpoints;
    damping = Determine_damping();
    
    cpoints = new ControlPoint[numOfpoints];
    springs = new Spring[numOfsprings];
    
    for (int i = 0; i < numOfpoints; i++) cpoints[i] = new ControlPoint(this.coords.get(i), pointMass, window);
    for (int i = 0; i < numOfsprings; i++) springs[i] = new Spring( cpoints[i], cpoints[i+1], segLength, stiffness, damping, window);
    
    xcurrent = new float[numOfpoints]; xnew = new float[numOfpoints];
    ycurrent = new float[numOfpoints]; ynew = new float[numOfpoints];
    vxcurrent = new float[numOfpoints]; vxnew = new float[numOfpoints];
    vycurrent = new float[numOfpoints]; vynew = new float[numOfpoints];
    
    dtmax = Determine_time_step();
    
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
  
  // calculate the pressure force at each point
  //PVector pressForce ( Field p ) {
  //  int orthSize = orth.length;
    
  //  PVector [] pf = new PVector[numOfpoints];
    
  //  pf[0] = 
    
  //  for ( OrthoNormal o: orth ) {
  //    float pdl = p.linear( o.cen.x, o.cen.y )*o.l;
  //    pv.add(pdl*o.nx, pdl*o.ny, 0);
  //  }
  //  return pv;
  //}
  
  // Calculate damping coefficient for numerical purposes
  float Determine_damping() {
    float d = sqrt(stiffness*pointMass);
    return d;
  }
  
  // Determine the size of the time step
  float Determine_time_step() {
    float ReLam, ImLam, dt;
    float n = float(numOfpoints);
    
    float RootDet = 1-stiffness*pointMass*sq(Length)/(2*sq(damping)*sq(n));
    float fact = 4*damping*sq(n)/(pointMass*sq(Length));
    
    if (RootDet<0) {
      ReLam = -fact;
      ImLam = -fact*sqrt(-RootDet);
    }
    else {
      ReLam = fact*(-1-sqrt(RootDet));
      ImLam = 0.0;
    }
    
    dt = -2*ReLam/(sq(ReLam)+sq(ImLam));
    println("dt is :"+dt);
    return dt;
  }
  
  // Move for testing purposes
  void move() {
    for (ControlPoint cp : cpoints) {
      cp.velocity = new PVector(0.1, .2);
      cp.position.add( cp.velocity );
    }
    getOrth();
  }
  
  // Apply internal Forces to particles
  void ApplyAllForces() {
    for (ControlPoint p : cpoints) p.clearForce();
    for (Spring s : springs) s.applyAllForces();
  }
  
  // Update the state of the sheet 
  void UpdateState(float [] Xnew, float [] Ynew, float [] VXnew, float [] VYnew) {
    for (int i = 0; i < numOfpoints; i++) {
      cpoints[i].position = new PVector(Xnew[i], Ynew[i]);
      cpoints[i].velocity = new PVector(VXnew[i], VYnew[i]);
    }
    getOrth();
  }
  
  // Trapezoid (Predictor-Corrector) Scheme
  void update(float dt, PVector ExtForce) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    
    int N = numOfpoints;
    
    float pMass = pointMass;
    PVector [] accelCurrent = new PVector[N];
    PVector [] accelPred = new PVector[N];
    
    float [] xPred = new float[N];
    float [] yPred = new float[N];
    float [] vxPred = new float[N];
    float [] vyPred = new float[N];
    
  //  // Store old state of sheet
  //  for (int i = 0; i < N; i++) {
  //    OldPosition[i] = prtcl[i].position.copy();
  //    OldVelocity[i] = prtcl[i].velocity.copy();
  //    prtcl[i].updatePositionOLD();
  //    prtcl[i].updateVelocityOLD();
  //  }
    
    // Apply Forces for this step
    ApplyAllForces();
    
    // calculate acceleration
    for (int i = 0; i < N; i++) {
      accelCurrent[i] = cpoints[i].force.copy();
      accelCurrent[i].div(pMass);
    }
    // accumulate any acceleration due to external forces
    for (int i = 0; i < N; i++) {
      accelCurrent[i].add(ExtForce);
    }
    
    // Calculate estimation
    for (int i = 0; i < N; i++) {
      if (!cpoints[i].fixed) {
        xPred[i] = xcurrent[i] + dt*vxcurrent[i];
        yPred[i] = ycurrent[i] + dt*vycurrent[i];
        vxPred[i] = vxcurrent[i] + dt*accelCurrent[i].x;
        vyPred[i] = vycurrent[i] + dt*accelCurrent[i].y;
      }
      else {
        xPred[i] = xcurrent[i];
        yPred[i] = ycurrent[i];
        vxPred[i] = vxcurrent[i];
        vyPred[i] = vycurrent[i];
      }
    }
    // Update the state of the filament for the correction
    UpdateState(xPred, yPred, vxPred, vyPred);
    
    // Apply Forces for the correction
    ApplyAllForces();
    
    // calculate acceleration for the correction
    for (int i = 0; i < N; i++) {
      accelPred[i] = cpoints[i].force.copy();
      accelPred[i].div(pMass);
    }
    // accumulate any acceleration due to external forces
    for (int i = 0; i < N; i++) {
      accelPred[i].add(ExtForce);
    }
    
    // Calculate at the new state
    for (int i = 0; i < N; i++) {
      if (!cpoints[i].fixed) {
        xnew[i] = xcurrent[i] + 0.5*dt*(vxcurrent[i]+vxPred[i]);
        ynew[i] = ycurrent[i] + 0.5*dt*(vycurrent[i]+vyPred[i]);
        vxnew[i] = vxcurrent[i] + 0.5*dt*(accelCurrent[i].x+accelPred[i].x);
        vynew[i] = vycurrent[i] + 0.5*dt*(accelCurrent[i].y+accelPred[i].y);
      }
      else {
        xnew[i] = xcurrent[i];
        ynew[i] = ycurrent[i];
        vxnew[i] = vxcurrent[i];
        vynew[i] = vycurrent[i];
      }
    }
    // Update the state of the filament for the correction
    UpdateState(xnew, ynew, vxnew, vynew);
  } // end of Trapezoid

} //=========== end of FlexibleSheet class ===============