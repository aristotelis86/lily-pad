//==================== FlexibleSheet Class ==================//



class FlexibleSheet extends LineSegBody {
  //============================= Attributes =================================//
  int numOfpoints;
  int numOfsprings;
  float Length, Mass, segLength, pointMass, stiffness, damping;
  float [] xcurrent, ycurrent, vxcurrent, vycurrent;
  PVector [] accelCurrent, accelPred;
  float [] xPred, yPred, vxPred, vyPred;
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
    
    for (int i = 0; i < numOfpoints; i++) cpoints[i] = new ControlPoint(this.coords.get(i), pointMass, thk, window);
    for (int i = 0; i < numOfsprings; i++) springs[i] = new Spring( cpoints[i], cpoints[i+1], segLength, stiffness, damping, window);
    
    dtmax = Determine_time_step();
    InitUpdateVars();
    
  } // end of Constructor
  
  //============================= Methods =================================//
  void InitUpdateVars() {
    xcurrent = new float[numOfpoints]; xPred = new float[numOfpoints]; xnew = new float[numOfpoints];
    ycurrent = new float[numOfpoints]; yPred = new float[numOfpoints]; ynew = new float[numOfpoints];
    vxcurrent = new float[numOfpoints]; vxPred = new float[numOfpoints]; vxnew = new float[numOfpoints];
    vycurrent = new float[numOfpoints]; vyPred = new float[numOfpoints]; vynew = new float[numOfpoints];
    accelCurrent = new PVector[numOfpoints];
    accelPred = new PVector[numOfpoints];
  }
  
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
  PVector [] pressForcePoints ( Field p ) {
    
    int orthSize = orth.length;
    PVector [] pf = new PVector[numOfpoints];
    for (int i=0; i<numOfpoints; i++) pf[i] = new PVector(0, 0);
    
    for ( int s=-1; s<=1; s+=2 ) {
      
      for ( int j=0; j<orthSize; j++ ) {
        float pdl = p.linear( cpoints[j].position.x+0.5*s*thk*orth[j].nx, cpoints[j].position.y+0.5*s*thk*orth[j].ny )*orth[j].l;
        PVector pTemp = new PVector(s*pdl*orth[j].nx, s*pdl*orth[j].ny);
        pf[j].sub(pTemp);
        pf[j+1].sub(pTemp);
      }
    }
    for (int j=1; j<numOfpoints-1; j++) pf[j].div(2);
    return pf;
  }
  
  // calculate the pressure force at each point
  PVector [] stressForcePoints ( VectorField Vel, float nu ) {
    
    Field u = Vel.x;
    Field v = Vel.y;
    Field grad = new Field( Vel.n, Vel.m);
    
    int orthSize = orth.length;
    PVector [] ps = new PVector[numOfpoints];
    for (int i=0; i<numOfpoints; i++) ps[i] = new PVector(0, 0);
    
    for ( int s=-1; s<=1; s+=2 ) {
      for ( int j=0; j<orthSize; j++ ) {
        
      }
    }
    
    return ps;
  }
  
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
      cpoints[i].UpdatePosition(Xnew[i], Ynew[i]);
      cpoints[i].UpdateVelocity(VXnew[i], VYnew[i]);
    }
    getOrth();
    getBox();
  }
  // Update the state of the sheet 
  void UpdateState(float [] Xnew, float [] Ynew) {
    for (int i = 0; i < numOfpoints; i++) {
      cpoints[i].UpdatePosition(Xnew[i], Ynew[i]);
    }
    getOrth();
    getBox();
  }
  
  // Get current positions and velocities
  void getState() {
    for(int i=0; i<numOfpoints; i++) {
      xcurrent[i] = cpoints[i].position.x;
      ycurrent[i] = cpoints[i].position.y;
      vxcurrent[i] = cpoints[i].velocity.x;
      vycurrent[i] = cpoints[i].velocity.y;
    }
  }
  
  
  // Trapezoid (Predictor-Corrector) Scheme
  void update(float dt, Field p) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    
    int N = numOfpoints;
    PVector [] ExtForce = new PVector[N];
    float pMass = pointMass;
    
    getState();
    ExtForce = pressForcePoints ( p );
  
    // Apply Forces for this step
    ApplyAllForces();
    
    // calculate acceleration
    for (int i = 0; i < N; i++) {
      accelCurrent[i] = cpoints[i].force.copy();
      accelCurrent[i].div(pMass);
    }
    // accumulate any acceleration due to external forces
    for (int i = 0; i < N; i++) {
      PVector ext = ExtForce[i].copy();
      ext.div(pMass);
      accelCurrent[i].add(ext);
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
  } // end of update (prediction)
  
  
  void update2(float dt, Field p) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    
    int N = numOfpoints;
    PVector [] ExtForce = new PVector[N];
    float pMass = pointMass;
    
    ExtForce = pressForcePoints ( p );
    
    // Apply Forces for the correction
    ApplyAllForces();
    
    // calculate acceleration for the correction
    for (int i = 0; i < N; i++) {
      accelPred[i] = cpoints[i].force.copy();
      accelPred[i].div(pMass);
    }
    // accumulate any acceleration due to external forces
    for (int i = 0; i < N; i++) {
      PVector ext = ExtForce[i].copy();
      ext.div(pMass);
      accelPred[i].add(ext);
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