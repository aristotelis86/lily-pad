////==================== CollisionSolver Class ==================//

////*** Depends on Particle, Spring and FlexibleSheet classes ***//

//class CollisionSolver {
  
//  //================ Attributes ================//
//  int N, Ns, NFil;
//  int [] NMass, NSpring;
//  FlexibleSheet FSheet;
  
//  ArrayList<ControlPoint> ListMass = new ArrayList<ControlPoint>();
//  ArrayList<Spring> ListSpring = new ArrayList<Spring>();
  
//  ControlPoint [] LocalMass;
//  Spring [] LocalSpring;
  
//  ArrayList<ControlPoint> NBoundPointCol; // For point-boundary collisions
//  ArrayList<ControlPoint> SBoundPointCol; // For point-boundary collisions
//  ArrayList<ControlPoint> WBoundPointCol; // For point-boundary collisions
//  ArrayList<ControlPoint> EBoundPointCol; // For point-boundary collisions
  
//  ArrayList<ControlPoint> PointACol; // For point-point collisions
//  ArrayList<ControlPoint> PointBCol; // For point-point collisions
  
//  ArrayList<ControlPoint> PointCol; // For point-edge collisions
//  ArrayList<Spring> SpringCol;  // For point-edge collisions
//  FloatList TimeCol;            // For point-edge collisions
//  FloatList LenCol;             // For point-edge collisions
//  float [] TR = new float[2];   // For point-edge collisions store time
//  float [] LR = new float[2];   // For point-edge collisions store location
//  int idx; // For point-edge collisions store index 
  
  
//  //================ Constructor ================//
//  CollisionSolver(FlexibleSheet filament_) { // single Sheet
    
//    FSheet = filament_;
//    N = FSheet.numOfpoints;
//    Ns = FSheet.numOfsprings;
    
//    LocalMass = new ControlPoint[N];
//    LocalSpring = new Spring[Ns];
    
//    LocalMass = FSheet.cpoints;
//    LocalSpring = FSheet.springs;
    
//    NBoundPointCol = new ArrayList<ControlPoint>();
//    SBoundPointCol = new ArrayList<ControlPoint>();
//    WBoundPointCol = new ArrayList<ControlPoint>();
//    EBoundPointCol = new ArrayList<ControlPoint>();
    
//    PointACol = new ArrayList<ControlPoint>();
//    PointBCol = new ArrayList<ControlPoint>();
    
//    PointCol = new ArrayList<ControlPoint>();
//    SpringCol = new ArrayList<Spring>();
//    TimeCol = new FloatList();
//    LenCol = new FloatList();
        
//  } // end of constructor #1
  
//  CollisionSolver(FlexibleSheet [] filament_) { // multiple sheets
//    NFil = filament_.length;
    
//    NMass = new int[NFil];
//    NSpring = new int[NFil];
    
//    for (int j = 0; j < NFil; j++) {
//      NMass[j] = filament_[j].numOfpoints;
//      NSpring[j] = filament_[j].numOfsprings;
//    }
    
//    N = 0; 
//    Ns = 0;
//    for (int nm : NMass) N += nm;
//    for (int ns : NSpring) Ns += ns;
    
//    for (int j = 0; j < NFil; j++) {
//      for (int i = 0; i < NMass[j]; i++) {
//        ListMass.add(filament_[j].cpoints[i]);
//      }
//    }
//    for (int j = 0; j < NFil; j++) {
//      for (int i = 0; i < NSpring[j]; i++) {
//        ListSpring.add(filament_[j].springs[i]);
//      }
//    }
    
//    LocalMass = new ControlPoint[N];
//    LocalSpring = new Spring[Ns];
    
//    for (int i = 0; i < N; i++) LocalMass[i] = ListMass.get(i);
//    for (int i = 0; i < Ns; i++) LocalSpring[i] = ListSpring.get(i);
    
//    NBoundPointCol = new ArrayList<ControlPoint>();
//    SBoundPointCol = new ArrayList<ControlPoint>();
//    WBoundPointCol = new ArrayList<ControlPoint>();
//    EBoundPointCol = new ArrayList<ControlPoint>();
    
//    PointACol = new ArrayList<ControlPoint>();
//    PointBCol = new ArrayList<ControlPoint>();
    
//    PointCol = new ArrayList<ControlPoint>();
//    SpringCol = new ArrayList<Spring>();
//    TimeCol = new FloatList();
//    LenCol = new FloatList();
    
//  } // end of constructor #2
  
  
//  //======================== Methods ===========================//
  
//  // Master Method - General Handling
//  void SolveCollisions() {
//    int collisionFlag = 1;
//    int pbFlag; // point-boundary collision counter
//    int ppFlag; // point-point collision counter
//    int peFlag; // point-edge collision counter
    
//    int iter = 0;
    
//    while ((boolean(collisionFlag)) && (iter<100)) {
//      iter++;
//      pbFlag = DetectBoundaryCollision();
//      ppFlag = DetectPointPointCollision();
//      peFlag = DetectPointEdgeCollision();
      
//      collisionFlag = pbFlag + ppFlag + peFlag;
      
//      if (boolean(collisionFlag)) {
//        // Resolve Collisions
//        ResolveBoundary();
//        ResolvePointPoint();
//        ResolvePointEdge();
//        //collisionFlag = 0;
//      }
//    }
//  } // end of SolveCollisions method
  
  
//  // Detect Boundary Collisions
//  int DetectBoundaryCollision() {
//    float clearRad;
//    int col_count = 0;
    
//    for (int i=0; i < N; i++) {
//      ControlPoint P = LocalMass[i];
//      clearRad = P.thick/2;
      
//      if (P.position.x < clearRad) {
//        col_count += 1;
//        WBoundPointCol.add(P);
//      }
//      else if (P.position.x > width - clearRad) {
//        col_count += 1;
//        EBoundPointCol.add(P);
//      }
//      if (P.position.y < clearRad) {
//        col_count += 1;
//        NBoundPointCol.add(P);
//      }
//      else if  (P.position.y > height - clearRad) {
//        col_count += 1;
//        SBoundPointCol.add(P);
//      }
//    } // end for loop over particles
    
//    return col_count;
//  } // end of Detect Boundary Collisions
  
//  // Detect Particle-Particle Collisions
//  int DetectPointPointCollision() {
//    float clearRad;
//    int col_count = 0;
    
//    for (int i = 0; i < N-1; i++) {
//      ControlPoint pi = LocalMass[i];
      
//      for (int j = i+1; j < N; j++) {
//        ControlPoint pj = LocalMass[j];
        
//        clearRad = (pi.thick + pj.thick)/6;
        
//        if (pi.position.dist(pj.position)<=clearRad) {
//          col_count += 1;
//          PointACol.add(pi);
//          PointBCol.add(pj);
//        }
//      } // end for loop over particles #2
//    } // end for loop over particles #1
    
//    return col_count;
    
//  } // end of Detect Point-Point Collisions
  
  
//  // Detect Particle-Spring Collisions
//  int DetectPointEdgeCollision() {
//    int col_count = 0;
    
//    for (int i = 0; i<N; i++) {
//      ControlPoint p = LocalMass[i];
//      for (int j=0; j<Ns; j++) {
//        Spring s = LocalSpring[j];
        
//        // continue only if the particle is not connected to the current spring
//        if ((p!=s.p1) && (p!=s.p2)) {
//          PVector C0 = p.positionOld.copy();
//          PVector C1 = p.position.copy();
          
//          PVector P0 = s.p1.positionOld.copy();
//          PVector P1 = s.p1.position.copy();
          
//          PVector Q0 = s.p2.positionOld.copy();
//          PVector Q1 = s.p2.position.copy();
          
//          PVector C0P0 = PVector.sub(C0,P0);
//          PVector C0Q0 = PVector.sub(C0,Q0);
          
//          PVector C1P1 = PVector.sub(C1,P1);
//          PVector C1Q1 = PVector.sub(C1,Q1);
          
//          PVector crossOld = C0P0.cross(C0Q0);
//          PVector crossNew = C1P1.cross(C1Q1);
          
//          if (crossOld.z*crossNew.z<0) {
//            boolean col;
//            col = Point_Edge_Penetration(P0, P1, Q0, Q1, C0, C1);
//            if (col) {
//              col_count += 1;
//              PointCol.add(p);
//              SpringCol.add(s);
//              TimeCol.append(TR[idx]);
//              LenCol.append(LR[idx]);
//            }
//          }
//        } // end if check for particle belonging to the spring
//      } // end of loop over springs
//    } // end of loop over particles
    
//    return col_count;
//  } // end of Detect Point-Edge Collisions
  
//  // Detect Point penetrating line segment
//  boolean Point_Edge_Penetration(PVector p0_, PVector p1_, PVector q0_, PVector q1_, PVector c0_, PVector c1_) {
    
//    boolean penetFlag = false;
//    float t, s;
    
//    PVector C0P0 = PVector.sub(c0_,p0_);
//    PVector Q0P0 = PVector.sub(q0_,p0_);
    
//    PVector C1C0P1P0 = PVector.sub(PVector.add(c1_,p0_),PVector.add(c0_,p1_));
//    PVector Q1Q0P1P0 = PVector.sub(PVector.add(q1_,p0_),PVector.add(q0_,p1_));
    
//    PVector A0 = C1C0P1P0.cross(Q1Q0P1P0);
//    PVector A11 = C1C0P1P0.cross(Q0P0);
//    PVector A12 = C0P0.cross(Q1Q0P1P0);
//    PVector A1 = PVector.add(A11,A12);
//    PVector A2 = C1C0P1P0.cross(Q1Q0P1P0);
    
//    TR = SolveQuadratic(A2.z,A1.z,A0.z);
    
//    for (int ii=0; ii<2; ii++) {
//      t = TR[ii];
      
//      PVector pt = PVector.add(p0_,PVector.mult(PVector.sub(p1_,p0_),t));
//      PVector qt = PVector.add(q0_,PVector.mult(PVector.sub(q1_,q0_),t));
//      PVector ct = PVector.add(c0_,PVector.mult(PVector.sub(c1_,c0_),t));
      
//      s = PVector.dot(PVector.sub(ct,pt),PVector.sub(qt,pt))/PVector.dot(PVector.sub(qt,pt),PVector.sub(qt,pt));
      
//      LR[ii] = s;
//    } // end for loop over time step ratios
    
//    if (((TR[0]>=0) && (LR[0]>=0)) && ((TR[0]<=1) && (LR[0]<=1))) {
//      penetFlag = true;
//      idx = 0;
//    }
//    else if (((TR[1]>=0) && (LR[1]>=0)) && ((TR[1]<=1) && (LR[1]<=1))) {
//      penetFlag = true;
//      idx = 1;
//    }
//    else penetFlag = false;
    
//    return penetFlag;
//  } // end of point-edge penetration detection method
  
  
//  // Solve Quadratic equation given the coefficients
//  float [] SolveQuadratic(float a, float b, float c) {
    
//    float [] Tf = new float[2]; 
    
//    float d = sq(b) - 4*a*c;
    
//    if ((a==0) && (b!=0)) {
//      Tf[0] = -c/b;
//      Tf[1] = -c/b;
//    }
//    else if ((a==0) && (b==0)) {
//      Tf[0] = 0.0;
//      Tf[1] = 0.0;
//    }
//    else {
//      if (d<0) {
//        Tf[0] = -999;
//        Tf[1] = -999;
//      }
//      else if (d==0) {
//        Tf[0] = -b/(2*a);
//        Tf[1] = -b/(2*a);
//      }
//      else {
//        Tf[0] = (-b+sqrt(d))/(2*a);
//        Tf[1] = (-b-sqrt(d))/(2*a);
//      }
//    }
    
//    if (Tf[0]>Tf[1]) {
//      float temp = Tf[0];
//      Tf[0] = Tf[1];
//      Tf[1] = temp;
//    }
    
//    return Tf;
//  } // end of quadratic solver
  
  
  
//  // Resolve Collisions with Boundaries
//  void ResolveBoundary() {
//    int NBP = NBoundPointCol.size();
//    int SBP = SBoundPointCol.size();
//    int WBP = WBoundPointCol.size();
//    int EBP = EBoundPointCol.size();
    
//    if (NBP>0) {
//      // Resolve north bound
//      for (int j = 0; j<NBP; j++) {
//        ControlPoint myP = NBoundPointCol.get(j);
//        myP.UpdatePosition(myP.position.x, myP.thick);
//        myP.UpdateVelocity(myP.velocity.x, -myP.velocity.y);
//      }
//      NBoundPointCol = new ArrayList<ControlPoint>();
//    }
//    if (SBP>0) {
//      // Resolve south bound
//      for (int j = 0; j<SBP; j++) {
//        ControlPoint myP = SBoundPointCol.get(j);
//        myP.UpdatePosition(myP.position.x, height-myP.thick);
//        myP.UpdateVelocity(myP.velocity.x, -myP.velocity.y);
//      }
//      SBoundPointCol = new ArrayList<ControlPoint>();
//    }
//    if (WBP>0) {
//      // Resolve west bound
//      for (int j = 0; j<WBP; j++) {
//        ControlPoint myP = WBoundPointCol.get(j);
//        myP.UpdatePosition(myP.thick, myP.position.y);
//        myP.UpdateVelocity(-myP.velocity.x, myP.velocity.y);
//      }
//      WBoundPointCol = new ArrayList<ControlPoint>();
//    }
//    if (EBP>0) {
//      // Resolve east bound
//      for (int j = 0; j<EBP; j++) {
//        ControlPoint myP = EBoundPointCol.get(j);
//        myP.UpdatePosition(width-myP.thick, myP.position.y);
//        myP.UpdateVelocity(-myP.velocity.x, myP.velocity.y);
//      }
//      EBoundPointCol = new ArrayList<ControlPoint>();
//    }
    
//  } // end of ResolveBoundary method
  
  
//  // Resolve Collisions occuring between points
//  void ResolvePointPoint() {
//    int Npp = PointACol.size();
//    float dtRatio = 0.5;
    
//    if (Npp>0) {
//      // resolve point-point
//      for (int j = 0; j<Npp; j++) {
//        ControlPoint pi = PointACol.get(j);
//        ControlPoint pj = PointBCol.get(j);
//        float piMass = pi.mass;
//        float pjMass = pj.mass;
        
//        float Vxi = (pi.velocity.x*(piMass-pjMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*pj.velocity.x;
//        float Vyi = (pi.velocity.y*(piMass-pjMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*pj.velocity.y;
//        float xinew = pi.positionOld.x + dtRatio*(pi.position.x-pi.positionOld.x);
//        float yinew = pi.positionOld.y + dtRatio*(pi.position.y-pi.positionOld.y);

//        float Vxj = (pj.velocity.x*(pjMass-piMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*pi.velocity.x;
//        float Vyj = (pj.velocity.y*(pjMass-piMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*pi.velocity.y;
//        float xjnew = pj.positionOld.x + dtRatio*(pj.position.x-pj.positionOld.x);
//        float yjnew = pj.positionOld.y + dtRatio*(pj.position.y-pj.positionOld.y);
        
//        if ((!pi.fixed) && (!pj.fixed)) {
//          pi.UpdateVelocity(Vxi, Vyi);
//          pj.UpdateVelocity(Vxj, Vyj);
//          pi.UpdatePosition(xinew,yinew);
//          pj.UpdatePosition(xjnew,yjnew);
//        }
//        else if ((pi.fixed) && (!pj.fixed)) {
//          pj.UpdateVelocity((-1)*pj.velocity.x,(-1)*pj.velocity.y);
//          pj.UpdatePosition(xjnew,yjnew);
//        }
//        else if ((!pi.fixed) && (pj.fixed)) {
//          pi.UpdateVelocity((-1)*pi.velocity.x,(-1)*pi.velocity.y);
//          pi.UpdatePosition(xinew,yinew);
//        }
//      }
//      PointACol = new ArrayList<ControlPoint>();
//      PointBCol = new ArrayList<ControlPoint>();
//    }
    
//  } // end ResolvePointPoint method
  
  
//  // Resolve Collisions between points and edges
//  void ResolvePointEdge() {
//    int Npe = PointCol.size();
//    float t, l, dtRatio;
//    float tol = 0.1;
    
//    if (Npe>0) {
//      // resolve point-edge
//      for (int j = 0; j<Npe; j++) {
//        ControlPoint p = PointCol.get(j);
//        Spring s = SpringCol.get(j);
//        t = TimeCol.get(j);
//        l = LenCol.get(j);
//        dtRatio = t-tol*t;
        
//        if (l>=0.5) {
//          // update s.p2 and p
//          float piMass = p.mass;
//          float pjMass = s.p2.mass;
          
//          float Vxi = (p.velocity.x*(piMass-pjMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*s.p2.velocity.x;
//          float Vyi = (p.velocity.y*(piMass-pjMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*s.p2.velocity.y;
//          float xinew = p.positionOld.x + dtRatio*(p.position.x-p.positionOld.x);
//          float yinew = p.positionOld.y + dtRatio*(p.position.y-p.positionOld.y);
  
//          float Vxj = (s.p2.velocity.x*(pjMass-piMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*p.velocity.x;
//          float Vyj = (s.p2.velocity.y*(pjMass-piMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*p.velocity.y;
//          float xjnew = s.p2.positionOld.x + dtRatio*(s.p2.position.x-s.p2.positionOld.x);
//          float yjnew = s.p2.positionOld.y + dtRatio*(s.p2.position.y-s.p2.positionOld.y);
          
//          if ((!p.fixed) && (!s.p2.fixed)) {
//            p.UpdateVelocity(Vxi, Vyi);
//            s.p2.UpdateVelocity(Vxj, Vyj);
//            p.UpdatePosition(xinew,yinew);
//            s.p2.UpdatePosition(xjnew,yjnew);
//          }
//          else if ((p.fixed) && (!s.p2.fixed)) {
//            s.p2.UpdateVelocity((-1)*s.p2.velocity.x,(-1)*s.p2.velocity.y);
//            s.p2.UpdatePosition(xjnew,yjnew);
//          }
//          else if ((!p.fixed) && (s.p2.fixed)) {
//            p.UpdateVelocity((-1)*p.velocity.x,(-1)*p.velocity.y);
//            p.UpdatePosition(xinew,yinew);
//          }
//        }
//        else if (l<0.5) {
//          // update s.p1 and p
//          float piMass = p.mass;
//          float pjMass = s.p1.mass;
          
//          float Vxi = (p.velocity.x*(piMass-pjMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*s.p1.velocity.x;
//          float Vyi = (p.velocity.y*(piMass-pjMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*s.p1.velocity.y;
//          float xinew = p.positionOld.x + dtRatio*(p.position.x-p.positionOld.x);
//          float yinew = p.positionOld.y + dtRatio*(p.position.y-p.positionOld.y);
  
//          float Vxj = (s.p1.velocity.x*(pjMass-piMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*p.velocity.x;
//          float Vyj = (s.p1.velocity.y*(pjMass-piMass)/(piMass+pjMass)) + (2*pjMass/(piMass+pjMass))*p.velocity.y;
//          float xjnew = s.p1.positionOld.x + dtRatio*(s.p1.position.x-s.p1.positionOld.x);
//          float yjnew = s.p1.positionOld.y + dtRatio*(s.p1.position.y-s.p1.positionOld.y);
          
//          if ((!p.fixed) && (!s.p1.fixed)) {
//            p.UpdateVelocity(Vxi, Vyi);
//            s.p1.UpdateVelocity(Vxj, Vyj);
//            p.UpdatePosition(xinew,yinew);
//            s.p1.UpdatePosition(xjnew,yjnew);
//          }
//          else if ((p.fixed) && (!s.p1.fixed)) {
//            s.p1.UpdateVelocity((-1)*s.p1.velocity.x,(-1)*s.p1.velocity.y);
//            s.p1.UpdatePosition(xjnew,yjnew);
//          }
//          else if ((!p.fixed) && (s.p1.fixed)) {
//            p.UpdateVelocity((-1)*p.velocity.x,(-1)*p.velocity.y);
//            p.UpdatePosition(xinew,yinew);
//          }
//        }
//      }
//      PointCol = new ArrayList<ControlPoint>();
//      SpringCol = new ArrayList<Spring>();
//      TimeCol = new FloatList();
//      LenCol = new FloatList();
//    }
    
    
//  } // end of ResolvePointEdge
  
  
//} // end of Collision solver class