/**********************************************************************
      FlexibleSheet class: Creates the flexible 1D object 
      consisting of control points and springs. Properties
      related to its behaviour are also created.
      Methods to handle and update the object can be found 
      here.
      
Example code:
Window view; // window in h-units
FloodPlot plot;

boolean saveimg = false;

int nx = (int)pow(2,6); // x-dir
int ny = (int)pow(2,6); // y-dir

float Re = 100;
float L = nx/4.;
float thick = 1;
int resol = 1;
float M = 10;
float stiff = 100;
PVector lpos = new PVector( nx/3., 0.55*ny );
PVector align = new PVector(1, 0);

float t=0;
float dt;

BDIM flow;
FlexibleSheet sheet;

void settings(){
    size(600, 600);
}

void setup() {
  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  
  sheet = new FlexibleSheet(L, thick, M, resol, stiff, lpos, align, view);
  sheet.cpoints[0].makeFixed(); // pinning leading point
  
  plot = new FloodPlot(view); // standard window
  
  dt = sheet.dtmax;
  
  flow = new BDIM(nx, ny, dt, sheet, 1/Re, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
} // end of setup


void draw() {
  
  sheet.waveOnSheet( t, sheet.Length/10., 2 );
  flow.update(sheet);
  flow.update2();
  
  plot.display(flow.p);
  sheet.display();
  
  if (saveimg) saveFrame("movie/frame_######.png");
  t += dt;
}
**********************************************************************/

class FlexibleSheet extends LineSegBody {
  //============================= Attributes =================================//
  int numOfpoints; // number of control points 
  int numOfsprings; // number of springs
  float Length, Mass, segLength, pointMass, stiffness, damping;
  float dtmax; // maximum time step to be used 
  
  // Dependancy from other objects
  ControlPoint [] cpoints; // stores the control points of the system
  Spring [] springs; // stores the connecting springs
  
  FlexibleSheet( float L_, float th_, float M_, int resol, float stiffness_, PVector lpos, PVector align, Window window, boolean damp ) {
    super(lpos.x, lpos.y, window); 
    thk = th_;
    weight = window.pdx(thk);
    
    align.normalize(); // alignment axis of the sheet
    PVector X0;
    for ( int i=0; i <= L_; i+=resol ) {
      X0 = PVector.add( lpos, PVector.mult(align, i));
      super.add(X0.x, X0.y);
    }
    super.end(false);
    
    Length = this.coords.get(0).dist(this.coords.get(this.coords.size()-1)); // find the length
    pointMass = M_; // linear mass distribution (Mass number)
    stiffness = stiffness_; // stiffness of each spring
    numOfpoints = this.coords.size();
    numOfsprings = numOfpoints - 1;
    segLength = resol; // resting length of each spring
    Mass = pointMass * numOfpoints; // total mass of the sheet
    
    // define if damping is needed or not
    if (damp) {
      damping = Determine_damping(); 
    }
    else damping = 0;
    
    cpoints = new ControlPoint[numOfpoints];
    springs = new Spring[numOfsprings];
    
    for (int i = 0; i < numOfpoints; i++) cpoints[i] = new ControlPoint( this.coords.get(i), pointMass, thk/2, window );
    for (int i = 0; i < numOfsprings; i++) springs[i] = new Spring( cpoints[i], cpoints[i+1], segLength, stiffness, damping, thk, window );
    
    // find the maximum time step size
    if (damp) dtmax = Determine_time_step();
    else dtmax = EigStiffMatrix();
    
  } // end of Constructor
  
  FlexibleSheet( float L_, float th_, float M_, int resol, float stiffness_, PVector lpos, PVector align, Window window ) {
    this( L_, th_, M_, resol, stiffness_, lpos, align, window, true);
  }
  
  //============================= Methods =================================//
  // !!!! override parent class method !!!!
  // Provides info to the fluid solver that it's not stationary object 
  boolean unsteady() {return true;}
  
  // !!!! override parent class method !!!!
  // Velocity of the fluid at some point from the surface of the object 
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
  } // end velocity()
  
  // Display control points for testing
  void pointsDisplay() {
    for (ControlPoint cp : cpoints) cp.display();
  }
  
  // Calculate current length of sheet (due to stretching or similar)
  float CurrentLength() {
    float newL = 0;
    for (int i=0; i<numOfpoints-1; i++) newL += this.coords.get(i).dist(this.coords.get(i+1));
    return newL;
  }
  
  // Update the state of the sheet 
  void UpdateState( float [] Xnew, float [] Ynew, float [] VXnew, float [] VYnew ) {
    for (int i = 0; i < numOfpoints; i++) {
      cpoints[i].UpdatePosition(Xnew[i], Ynew[i]);
      cpoints[i].UpdateVelocity(VXnew[i], VYnew[i]);
    }
    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  }
  
  // Update the state of the sheet 
  void UpdateState( float [] Xnew, float [] Ynew ) {
    for (int i = 0; i < numOfpoints; i++) {
      cpoints[i].UpdatePosition(Xnew[i], Ynew[i]);
    }
    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  }
  
  // Update the position of the center of the Body class (xc)
  void UpdateCenter() {
    super.xc.x = cpoints[0].position.x;
    super.xc.y = cpoints[0].position.y;
  }
  
  // Determine new steady state based on constant external forcing
  void Calculate_Stretched_Positions( PVector F ) {
    float g = F.mag();
    float ll;
    PVector newll;
    float [] xnew = new float[numOfpoints];
    float [] ynew = new float[numOfpoints];
    PVector FDir = F.copy();
    FDir.normalize();
    
    xnew[0] = cpoints[0].position.x; ynew[0] = cpoints[0].position.y;
    for (int i = 0; i < numOfsprings; i++) {
      ll = (float(numOfpoints) - float(i) - 1)*(pointMass*g/stiffness) + segLength;
      newll = PVector.add(cpoints[i].position, PVector.mult(FDir,ll));
      cpoints[i+1].UpdatePosition(newll.x, newll.y);
    }
    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  } // end of Calculate_Stretched_Positions
  
  // Move for testing purposes (translation only)
  void move() {
    for (ControlPoint cp : cpoints) {
      cp.velocity = new PVector(0.1, .2);
      cp.position.add( cp.velocity );
    }
    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  }
  
  // Apply internal Forces to control points
  void ApplyIntForces() {
    for (Spring s : springs) s.ApplyAllForces();
  }
  
  // Apply external forces to control points
  void ApplyExtForces( PVector [] F) {
    for (int i=0; i<numOfpoints; i++) cpoints[i].force.add(F[i]);
  }
  
  // Apply gravitational Forces to particles
  void ApplyGravity( PVector g) {
    for (ControlPoint cp : cpoints) {
      PVector extg = g.copy();
      extg.mult(cp.mass);
      cp.ApplyForce(extg);
    }
  }
  
  // Clear all forces acting in control points
  void ClearForces() {
    for (ControlPoint p : cpoints) p.clearForce();
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
  
  // calculate the shear stress force at each point
  PVector [] stressForcePoints ( VectorField Vel, float nu ) {
    
    Field omega = Vel.curl();
    
    int orthSize = orth.length;
    PVector [] ps = new PVector[numOfpoints];
    for (int i=0; i<numOfpoints; i++) ps[i] = new PVector(0, 0);
    
    for ( int s=-1; s<=1; s+=2 ) {
      for ( int j=0; j<orthSize; j++ ) {
        float omegaZVal = omega.linear( cpoints[j].position.x+0.5*s*thk*orth[j].nx, cpoints[j].position.y+0.5*s*thk*orth[j].ny );
        PVector pTemp = new PVector(s*omegaZVal*orth[j].ny, s*omegaZVal*orth[j].nx);
        ps[j].add(pTemp);
        ps[j+1].add(pTemp);
      }
    }
    for (int j=1; j<numOfpoints-1; j++) ps[j].div(2);
    for (int j=0; j<numOfpoints; j++) ps[j].mult(nu);
    
    return ps;
  }
  
  // Update based on Predictor-Corrector Scheme (gravity only)
  void update(float dt, PVector p) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    
    // Apply Forces for this step
    ClearForces(); // clear forces from control points 
    ApplyIntForces(); // apply spring-damper forces
    ApplyGravity( p ); // apply gravity
    
    // run the update through control points
    for (ControlPoint cp : cpoints) cp.update( dt );

    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  } // end of update (prediction)
  
  void update2(float dt, PVector p) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    
    // Apply Forces for the correction
    ClearForces(); // clear forces from control points 
    ApplyIntForces(); // apply spring-damper forces
    ApplyGravity( p ); // apply gravity
    
    // run the update through control points
    for (ControlPoint cp : cpoints) cp.update2( dt );

    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  } // end of Trapezoid
  
  // Alternative Update based on Predictor-Corrector Scheme (gravity only)
  void updateAlt(float dt, PVector p) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    
    // Apply Forces for this step
    ClearForces(); // clear forces from control points 
    ApplyIntForces(); // apply spring-damper forces
    ApplyGravity( p ); // apply gravity
    
    // run the update through control points
    for (ControlPoint cp : cpoints) cp.updateAlt( dt );

    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  } // end of update (alternative prediction, already 2nd order in space)
  
  void updateAlt2(float dt, PVector p) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    
    // Apply Forces for the correction
    ClearForces(); // clear forces from control points
    ApplyIntForces(); // apply spring-damper forces
    ApplyGravity( p ); // apply gravity
    
    // run the update through control points
    for (ControlPoint cp : cpoints) cp.updateAlt2( dt );
    
  } // end of alternative prediction-correction (2nd step update velocity only)
  
  /*************** Fluid present **********************************/
  // Update based on Predictor-Corrector Scheme
  void update( float dt, BDIM flow ) { update( dt, flow, new PVector(0,0)); }
  void update( float dt, BDIM flow, PVector g) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    PVector [] pressForce = new PVector[numOfpoints];
    PVector [] stressForce = new PVector[numOfpoints];
    
    pressForce = pressForcePoints ( flow.p ); // get pressure force
    stressForce = stressForcePoints( flow.u, flow.nu ); // get stress force
    
    // Apply Forces for this step
    ClearForces(); // clear forces from control points
    ApplyIntForces(); // apply spring-damper forces
    ApplyExtForces( pressForce ); // apply pressure
    ApplyExtForces( stressForce ); // apply stress
    ApplyGravity( g ); // apply gravity (if any)
    
    // run the update through control points
    for (ControlPoint cp : cpoints) cp.update( dt );

    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  } // end of update (prediction)
  
  void update2( float dt, BDIM flow ) { update2( dt, flow, new PVector(0,0)); }
  void update2(float dt, BDIM flow, PVector g) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    PVector [] pressForce = new PVector[numOfpoints];
    PVector [] stressForce = new PVector[numOfpoints];
    
    pressForce = pressForcePoints ( flow.p ); // get pressure force
    stressForce = stressForcePoints( flow.u, flow.nu ); // get stress force
    
    // Apply Forces for the correction
    ClearForces(); // clear forces from control points
    ApplyIntForces(); // apply spring-damper forces
    ApplyExtForces( pressForce ); // apply pressure
    ApplyExtForces( stressForce ); // apply stress
    ApplyGravity( g ); // apply gravity (if any)
    
    // run the update through control points
    for (ControlPoint cp : cpoints) cp.update2( dt );
      
    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  } // end of Trapezoid
  
  // Alternative Update based on Predictor-Corrector Scheme
  void updateAlt( float dt, BDIM flow ) { updateAlt( dt, flow, new PVector(0,0)); }
  void updateAlt( float dt, BDIM flow, PVector g) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    PVector [] pressForce = new PVector[numOfpoints];
    PVector [] stressForce = new PVector[numOfpoints];
    
    pressForce = pressForcePoints ( flow.p ); // get pressure force
    stressForce = stressForcePoints( flow.u, flow.nu ); // get stress force
    
    // Apply Forces for this step
    ClearForces(); // clear forces from control points
    ApplyIntForces(); // apply spring-damper forces
    ApplyExtForces( pressForce ); // apply pressure
    ApplyExtForces( stressForce ); // apply stress
    ApplyGravity( g ); // apply gravity (if any)
    
    // run the update through control points
    for (ControlPoint cp : cpoints) cp.updateAlt( dt );
    
    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  } // end of update (prediction)
  
  void updateAlt2( float dt, BDIM flow ) { updateAlt2( dt, flow, new PVector(0,0)); }
  void updateAlt2(float dt, BDIM flow, PVector g) {
    
    if (dt>dtmax) {
      println("WARNING dt constrained to maximum permitted:"+dtmax);
      dt = dtmax;
    }
    PVector [] pressForce = new PVector[numOfpoints];
    PVector [] stressForce = new PVector[numOfpoints];
    
    pressForce = pressForcePoints ( flow.p ); // get pressure force
    stressForce = stressForcePoints( flow.u, flow.nu ); // get stress
    
    // run the update through control points
    ClearForces(); // clear forces from control points
    ApplyIntForces(); // apply spring-damper forces
    ApplyExtForces( pressForce ); // apply pressure
    ApplyExtForces( stressForce ); // apply stress
    ApplyGravity( g ); // apply gravity (if any)
    
    // calculate acceleration for the correction
    for (ControlPoint cp : cpoints) cp.updateAlt2( dt );
  } // end of Trapezoid
  
  // Calculate damping coefficient for numerical purposes
  float Determine_damping() {
    float d = 2*sqrt(stiffness*pointMass);
    return d;
  }
  
  // Determine the size of the time step (damping is on)
  float Determine_time_step() {
    float ReLam, ImLam, dt;
    float n = float(numOfsprings);
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
    
    dt = -1*ReLam/(sq(ReLam)+sq(ImLam));
    //dt = -1/(sq(ReLam)+sq(ImLam));
    return dt;
  } // end of Determine_time_step
  
  
  // Find dominant eigenvalue of stiffness matrix
  float EigStiffMatrix() {
    float [][] K = new float[numOfpoints][numOfpoints];
    float [] X0 = new float[numOfpoints];
    float [] Xn;
    float lambda;
    
    for (int i=0; i<numOfpoints; i++) X0[i] = pow(-1,i);
    
    for (int i=0; i<numOfpoints; i++) {
      for (int j=0; j<numOfpoints; j++) {
        if (i==j) {
          if ((i==0) || (i==numOfpoints-1)) K[i][j] = -1;
          else K[i][j] = -2;
        }
        else if ((j==i+1) || (j==i-1)) K[i][j] = 1;
        else K[i][j] = 0; 
      }
    }
    
    Xn = MatProduct( K, X0);
    for (int n=0; n<10*numOfpoints; n++) {
      Xn = MatProduct( K, Xn);
    }
    
    lambda = ArrayProd( Xn, X0 )/ArrayProd( X0, X0);
    
    float dt = pointMass/abs(stiffness*lambda);
    return dt;
  }
  
  // Calculate the product of square matrices
  float [] MatProduct( float [][] A, float [] c ) {
    int iRow = A.length;
    int jColumn = A[0].length;
    int resSize = c.length;
    float [] result = new float[resSize];
    float [] absResult = new float[resSize];
    float minResult;
    for (int i=0; i<resSize; i++) result[i] = 0;
    
    for (int i=0; i<iRow; i++) {
      for (int j=0; j<jColumn; j++) {
        result[i] += A[i][j]*c[j];
      }
    }
    
    for (int i=0; i<resSize; i++) absResult[i] = abs(result[i]);
    minResult = absResult[0];
    for (int i=1; i<resSize; i++) {
      if (absResult[i]<minResult) minResult = absResult[i];
    }
    for (int i=0; i<resSize; i++) result[i] = result[i]/minResult;
    return result;
  }
  
  // Calculate product of arrays
  float ArrayProd( float [] a, float [] b ) {
    int resSize = a.length;
    float result = 0;
    for (int i=0; i<resSize; i++) result += a[i]*b[i]; 
    return result;
  }
  
  // Prescribed motion of the sheet to test effect on flow
  void waveOnSheet( float t, float sinAmp, float sinN ) {
    
    int nn = cpoints.length;
    float [] x = new float[nn];
    float [] y = new float[nn];
    float [] vx = new float[nn];
    float [] vy = new float[nn];
    
    x[0] = cpoints[0].position.x;
    y[0] = cpoints[0].position.y;
    vx[0] = 0.;
    vy[0] = 0.;
    for (int i = 1; i < nn; i++) {
      x[i] = (Length/(nn-1)) + x[i-1];
      y[i] = (sinAmp * sin(sinN*PI*(x[i]-x[0])/Length))*sin(2*t) + y[0];
      vx[i] = 0.;
      vy[i] = 2*cos(2*t)*(sinAmp * sin(sinN*PI*(x[i]-x[0])/Length));
    }
    UpdateState(x, y, vx, vy);
    this.UpdateCenter();
    super.getOrth();
    super.getBox();
  }
  
} //=========== end of FlexibleSheet class ===============
