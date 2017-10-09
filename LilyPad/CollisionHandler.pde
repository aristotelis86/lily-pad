/**********************************************************************
      CollisionHandler class: ControlPoints and/or FlexibleSheets
      are gathered by this class to handle any collisions during
      the simulation. Detection and resolution are separated.

**********************************************************************/
class CollisionHandler {
  //================= Attributes ====================//
  int Ncp; // number of control points to consider
  int Nsp; // number of springs to consider
  
  ArrayList<ControlPoint> LocalCPoints = new ArrayList<ControlPoint>(); // gather all control points
  ArrayList<Spring> LocalSprings = new ArrayList<Spring>(); // gather all springs
  
  //================= Constructor ====================//
  CollisionHandler() {
    Ncp = 0;
    Nsp = 0;
    // Just exploit the class for its methods
  }
  CollisionHandler( ControlPoint [] cpoints ) {
    Ncp = cpoints.length;
    Nsp = 0;
    for (int i=0; i<Ncp; i++) LocalCPoints.add(cpoints[i]);
  }
  CollisionHandler( ControlPoint cpoint ) {
    Ncp = 1;
    Nsp = 0;
    LocalCPoints.add(cpoint);
  }
  CollisionHandler( ControlPoint [] cpoints, Spring [] springs ) {
    Ncp = cpoints.length;
    Nsp = springs.length;
    for (int i=0; i<Ncp; i++) LocalCPoints.add(cpoints[i]);
    for (int i=0; i<Nsp; i++) LocalSprings.add(springs[i]);
  }
  CollisionHandler( ControlPoint [] cpoints, Spring spring ) {
    Ncp = cpoints.length;
    Nsp = 1;
    for (int i=0; i<Ncp; i++) LocalCPoints.add(cpoints[i]);
    LocalSprings.add(spring);
  }
  CollisionHandler( FlexibleSheet [] sheet ) {
    int Nfs = sheet.length;
    
    Ncp = 0;
    Nsp = 0;
    for (int ij=0; ij<Nfs; ij++) {
      for (int i=0; i<sheet[ij].numOfpoints; i++) {
        LocalCPoints.add(sheet[ij].cpoints[i]);
        Ncp++;
      }
      for (int i=0; i<sheet[ij].numOfsprings; i++) {
        LocalSprings.add(sheet[ij].springs[i]);
        Nsp++;
      }
    }
  }
  CollisionHandler( FlexibleSheet sheet ) {
    Ncp = sheet.numOfpoints;
    Nsp = sheet.numOfsprings;
    for (int i=0; i<Ncp; i++) {
      LocalCPoints.add(sheet.cpoints[i]);
    }
    for (int i=0; i<Nsp; i++) {
      LocalSprings.add(sheet.springs[i]);
    }
  }
  
  //================= Method ====================//
  // Detect Boundary Collisions
  void DetectBoundCollision() {
    for (int i=0; i<Ncp; i++) {
      ControlPoint cp = LocalCPoints.get(i);
      if (cp.position.x < cp.diameter/2) {
        this.ResolveBoundCollisions( "West", cp );
      }
      else if (cp.position.x > cp.myWindow.x.inE - cp.diameter/2) {
        this.ResolveBoundCollisions( "East", cp );
      }
      if (cp.position.y < cp.diameter/2) {
        this.ResolveBoundCollisions( "North", cp );
      }
      else if (cp.position.y > cp.myWindow.y.inE - cp.diameter/2) {
        this.ResolveBoundCollisions( "South", cp );
      }
    }
  }
  // Detect Boundary Collisions
  void DetectBoundCollision( ControlPoint cp ) {
    if (cp.position.x < cp.diameter/2) {
      this.ResolveBoundCollisions( "West", cp );
    }
    else if (cp.position.x > cp.myWindow.x.inE - cp.diameter/2) {
      this.ResolveBoundCollisions( "East", cp );
    }
    if (cp.position.y < cp.diameter/2) {
      this.ResolveBoundCollisions( "North", cp );
    }
    else if (cp.position.y > cp.myWindow.y.inE - cp.diameter/2) {
      this.ResolveBoundCollisions( "South", cp );
    }
  }
  // Resolve boundary collisions
  void ResolveBoundCollisions( String bound, ControlPoint cp ){
    float e = 0.95;
    if ((bound.equals("North")==true) || (bound.equals("South")==true)) {
      float vx = cp.velocity.x;
      float vy = -e*cp.velocity.y;
      float x = cp.position.x;
      float y = cp.position.y;
      if (bound.equals("North")==true) y = cp.diameter/2;
      else if (bound.equals("South")==true) y = cp.myWindow.y.inE-cp.diameter/2;
      cp.UpdatePosition(x,y); cp.UpdateVelocity(vx,vy);
    }
    if ((bound.equals("West")==true) || (bound.equals("East")==true)) {
      float vx = -e*cp.velocity.x;
      float vy = cp.velocity.y;
      float x = cp.position.x;
      float y = cp.position.y;
      if (bound.equals("West")==true) x = cp.diameter/2;
      else if (bound.equals("East")==true) x = cp.myWindow.x.inE-cp.diameter/2;
      cp.UpdatePosition(x,y); cp.UpdateVelocity(vx,vy);
    }
  }
  
  // Detect CPoint-CPoint collisions
  void DetectCPointCPointCollision() {
    int N = LocalCPoints.size();
    if (N>1) {
      for (int i=0; i<N-1; i++) {
        for (int j=i+1; j<N; j++) {
          ControlPoint cpi = LocalCPoints.get(i);
          ControlPoint cpj = LocalCPoints.get(j);
          float dij = cpi.distance(cpj);
          float clearRad = (cpi.diameter + cpj.diameter)/2.;
          if (dij<=clearRad) {
            float penet = 0.5*(cpi.diameter + cpj.diameter) - dij;
            float rewindi = penet/(0.5*cpi.diameter);
            float rewindj = penet/(0.5*cpj.diameter);
            this.ResolveCPCPCollisions( cpi, rewindi, cpj, rewindj );
            // Must include a detection method for fast moving cpoints
            // Cheking if their relative velocity is greater than clearRad 
          }
        }
      }
    }
  }
  // Detect CPoint-CPoint collisions
  boolean DetectCPointCPointCollision( ControlPoint cpi, ControlPoint cpj ) {
    boolean Flag = false;
    float dij = cpi.distance(cpj);
    float clearRad = (cpi.diameter + cpj.diameter)/2.;
    if (dij<=clearRad) {
      float penet = 0.5*(cpi.diameter + cpj.diameter) - dij;
      float rewindi = penet/(0.5*cpi.diameter);
      float rewindj = penet/(0.5*cpj.diameter);
      this.ResolveCPCPCollisions( cpi, rewindi, cpj, rewindj );
      Flag = true;
      // Must include a detection method for fast moving cpoints
     // Cheking if their relative velocity is greater than clearRad 
    }
    return Flag;
  }
  // Resolve CPoint-CPoint collisions
  void ResolveCPCPCollisions( ControlPoint cpi,  float rewti, ControlPoint cpj, float rewtj ) {
    cpi.rewindPosition(rewti);
    cpj.rewindPosition(rewtj);
    float Vxi = (cpi.velocity.x*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.x;
    float Vyi = (cpi.velocity.y*(cpi.mass-cpj.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpj.velocity.y;
    float Vxj = (cpj.velocity.x*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.x;
    float Vyj = (cpj.velocity.y*(cpj.mass-cpi.mass)/(cpi.mass+cpj.mass)) + (2*cpj.mass/(cpi.mass+cpj.mass))*cpi.velocity.y;
    cpi.UpdateVelocity( Vxi, Vyi );
    cpj.UpdateVelocity( Vxj, Vyj );
  }
  
  // Detect spring-spring intersection
  void DetectSpringSpringCollision() {
    if (Nsp>1) {
      for (int i=0; i<Nsp-1; i++) {
        for (int j=i+1; j<Nsp; j++) {
          Spring spi = LocalSprings.get(i);
          Spring spj = LocalSprings.get(j);
          ControlPoint cpi1 = spi.p1;
          ControlPoint cpj1 = spj.p1;
          ControlPoint cpi2 = spi.p2;
          ControlPoint cpj2 = spj.p2;
          if ((((cpi1!=cpj1) && (cpi1!=cpj2))) &&  (((cpi2!=cpj1) && (cpi2!=cpj2)))) {
            BoundBox Boxi, Boxj;
            Boxi = new BoundBox( spi );
            Boxj = new BoundBox( spj );
            //Boxi.displayOBB();
            //Boxj.displayOBB();
            boolean flag = this.doBoundBoxesOverlap( Boxi, Boxj );
            if (flag) {
              this.ResolveSpringSpringCollisions( spi, spj );
              //noLoop();
            }
          }
        }
      }
    }
  }
  void ResolveSpringSpringCollisions( Spring sp1, Spring sp2 ) {
    
      float sp1Velx = 0.5*(sp1.p1.velocity.x + sp1.p2.velocity.x);
      float sp1Vely = 0.5*(sp1.p1.velocity.y + sp1.p2.velocity.y);
      float sp1Mass = sp1.p1.mass + sp1.p2.mass;
      
      float sp2Velx = 0.5*(sp2.p1.velocity.x + sp2.p2.velocity.x);
      float sp2Vely = 0.5*(sp2.p1.velocity.y + sp2.p2.velocity.y);
      float sp2Mass = sp2.p1.mass + sp2.p2.mass;
      
      float Vx1 = (sp1Velx*(sp1Mass-sp2Mass)/(sp1Mass+sp2Mass)) + (2*sp2Mass/(sp1Mass+sp2Mass))*sp2Velx;
      float Vy1 = (sp1Vely*(sp1Mass-sp2Mass)/(sp1Mass+sp2Mass)) + (2*sp2Mass/(sp1Mass+sp2Mass))*sp2Vely;
      float Vx2 = (sp2Velx*(sp2Mass-sp1Mass)/(sp1Mass+sp2Mass)) + (2*sp2Mass/(sp1Mass+sp2Mass))*sp1Velx;
      float Vy2 = (sp2Vely*(sp2Mass-sp1Mass)/(sp1Mass+sp2Mass)) + (2*sp2Mass/(sp1Mass+sp2Mass))*sp1Vely;
      sp1.p1.UpdateVelocity( Vx1, Vy1 );
      sp1.p2.UpdateVelocity( Vx1, Vy1 );
      sp2.p1.UpdateVelocity( Vx2, Vy2 );
      sp2.p2.UpdateVelocity( Vx2, Vy2 );
  }
  
  boolean doBoundBoxesOverlap( BoundBox b1, BoundBox b2 ) {
    Boolean overFlag = false;
    PVector [] CheckAxes = new PVector[4];
    Boolean [] AxesColFlag = new Boolean[4];
    CheckAxes[0] = new PVector(b1.lines[0].orth.nx, b1.lines[0].orth.ny);
    CheckAxes[1] = new PVector(b1.lines[1].orth.nx, b1.lines[1].orth.ny);
    CheckAxes[2] = new PVector(b2.lines[0].orth.nx, b2.lines[0].orth.ny);
    CheckAxes[3] = new PVector(b2.lines[1].orth.nx, b2.lines[1].orth.ny);
    for (Boolean ax : AxesColFlag) {
      ax = false;
    }
    
    axesLoop : for (int i=0; i<4; i++) {
      float [] ProjB1 = new float[5];
      float [] ProjB2 = new float[5];
      for (int j=0; j<5; j++) {
        ProjB1[j] = PVector.dot(CheckAxes[i],b1.vertices[j]);
        ProjB2[j] = PVector.dot(CheckAxes[i],b2.vertices[j]);
      }
      float minProjB1 = ProjB1[1];
      float maxProjB1 = ProjB1[1];
      float minProjB2 = ProjB2[1];
      float maxProjB2 = ProjB2[1];
      
      for (int j=2; j<5; j++) {
        if (ProjB1[j]<minProjB1) minProjB1 = ProjB1[j];
        if (ProjB1[j]>maxProjB1) maxProjB1 = ProjB1[j];
        if (ProjB2[j]<minProjB2) minProjB2 = ProjB2[j];
        if (ProjB2[j]>maxProjB2) maxProjB2 = ProjB2[j];
      }
      if ((maxProjB2 < minProjB1) || (maxProjB1 < minProjB2)) {
        AxesColFlag[i] = false;
        break axesLoop;
      }
      else AxesColFlag[i] = true;
    }
    if (AxesColFlag[0] && AxesColFlag[1] && AxesColFlag[2] && AxesColFlag[3]) {
      overFlag = true;
    }
    return overFlag;
  }
  
  // Main Handling method
  void HandleCollisions() {
    // There should be a loop that restarts everytime a collision happened.
    // In this way a collision-free state is ensured every time...
    this.DetectBoundCollision();
    this.DetectCPointCPointCollision();
    this.DetectSpringSpringCollision();
  }
  
} // end of class
