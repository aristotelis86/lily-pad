class WriteInfo {

  PrintWriter [] outputPos; // output for positions
  PrintWriter [] outputVel; // output for velocities
  PrintWriter [] outputForce; // output for forces
  PrintWriter [] outputEnergy; // output for energy

  int Ncp, Nsg, Nfs; // number of lines to write per call
  ArrayList<ControlPoint> [] myCPoints;
  ArrayList<Spring> [] mySprings;

  boolean posFlag = false;
  boolean velFlag = false;
  boolean forFlag = false;
  boolean eneFlag = false;

  //================= Constructor ====================//
  WriteInfo(ControlPoint cp) {
    Ncp = 1;
    outputPos = new PrintWriter[1];
    outputVel = new PrintWriter[1];
    outputForce = new PrintWriter[1];
    outputEnergy = new PrintWriter[1];
    
    myCPoints =  new ArrayList[1];
    myCPoints[0] = new ArrayList<ControlPoint>();
    myCPoints[0].add(cp);
    InitPositionOut(); posFlag = true; InitPoints(outputPos[0]);
    InitVelocityOut(); velFlag = true; InitPoints(outputVel[0]);
    InitForceOut(); forFlag = true; InitPoints(outputForce[0]);
    InitEnergyOut(); eneFlag = true; InitPoints(outputEnergy[0]);
  }

  WriteInfo(ControlPoint [] cp) {
    Ncp = cp.length;
    outputPos = new PrintWriter[1];
    outputVel = new PrintWriter[1];
    outputForce = new PrintWriter[1];
    outputEnergy = new PrintWriter[1];
    
    myCPoints =  new ArrayList[1];
    myCPoints[0] = new ArrayList<ControlPoint>();
    
    for (int i=0; i<Ncp; i++) { myCPoints[0].add(cp[i]); }
    InitPositionOut(); posFlag = true; InitPoints(outputPos[0]);
    InitVelocityOut(); velFlag = true; InitPoints(outputVel[0]);
    InitForceOut(); forFlag = true; InitPoints(outputForce[0]);
    InitEnergyOut(); eneFlag = true; InitPoints(outputEnergy[0]);
  }

  WriteInfo( FlexibleSheet fs ) {
    Nfs = 1;
    Ncp = fs.numOfpoints;
    Nsg = fs.numOfsprings;
    
    outputPos = new PrintWriter[1];
    outputVel = new PrintWriter[1];
    outputForce = new PrintWriter[1];
    outputEnergy = new PrintWriter[1];
    
    myCPoints =  new ArrayList[1]; 
    mySprings =  new ArrayList[1];
    myCPoints[0] = new ArrayList<ControlPoint>();
    mySprings[0] = new ArrayList<Spring>();

    for (int i=0; i<Ncp; i++) { 
      myCPoints[0].add(fs.cpoints[i]);
    }
    for (int i=0; i<Nsg; i++) { 
      mySprings[0].add(fs.springs[i]);
    }
    
    InitPositionOut(); posFlag = true;
    InitVelocityOut(); velFlag = true;
    InitForceOut(); forFlag = true;
    InitEnergyOut(); eneFlag = true;

    InitSheetOut( fs );
  }

  WriteInfo(FlexibleSheet [] fs) {
    Nfs = fs.length;

    outputPos = new PrintWriter[Nfs];
    outputVel = new PrintWriter[Nfs];
    outputForce = new PrintWriter[Nfs];
    outputEnergy = new PrintWriter[Nfs];

    myCPoints =  new ArrayList[Nfs]; 
    mySprings =  new ArrayList[Nfs];


    for ( int j=0; j<Nfs; j++ ) {
      myCPoints[j] = new ArrayList<ControlPoint>();
      mySprings[j] = new ArrayList<Spring>();
      for (int i=0; i<fs[j].numOfpoints; i++) { 
        myCPoints[j].add(fs[j].cpoints[i]);
      }

      InitPositionOut( j ); posFlag = true;
      InitVelocityOut( j ); velFlag = true;
      InitForceOut( j ); forFlag = true;

      for (int i=0; i<fs[j].numOfsprings; i++) { 
        mySprings[j].add(fs[j].springs[i]);
      }

      InitEnergyOut( j ); eneFlag = true;
    }
    InitSheetOut( fs );
  }

  //================= Methods ====================//
  void InitPositionOut( int d) {
    outputPos[d] = createWriter("./info/positions"+d+".txt");
    outputPos[d].println("=========== Positions ==========");
  }
  void InitPositionOut() { 
    InitPositionOut( 0 );
  }

  void InitVelocityOut( int d) {
    outputVel[d] = createWriter("./info/velocities"+d+".txt");
    outputVel[d].println("=========== Velocities ==========");
  }
  void InitVelocityOut() { 
    InitVelocityOut( 0 );
  }

  void InitForceOut( int d ) {
    outputForce[d] = createWriter("./info/forces"+d+".txt");
    outputForce[d].println("=========== Forces ==========");
  }
  void InitForceOut() { 
    InitForceOut( 0 );
  }

  void InitEnergyOut( int d ) {
    outputEnergy[d] = createWriter("./info/energy"+d+".txt");
    outputEnergy[d].println("=========== Energy ==========");
  }
  void InitEnergyOut() { 
    InitEnergyOut( 0 );
  }

  void InitSheetOut( FlexibleSheet fs ) {
    outputPos[0].println("Length: "+fs.Length+" Mass: "+fs.Mass+" Points: "+fs.numOfpoints+" Stiffness: "+fs.stiffness+" Damping: "+fs.damping+" dt: "+fs.dtmax);
    outputVel[0].println("Length: "+fs.Length+" Mass: "+fs.Mass+" Points: "+fs.numOfpoints+" Stiffness: "+fs.stiffness+" Damping: "+fs.damping+" dt: "+fs.dtmax);
    outputForce[0].println("Length: "+fs.Length+" Mass: "+fs.Mass+" Points: "+fs.numOfpoints+" Stiffness: "+fs.stiffness+" Damping: "+fs.damping+" dt: "+fs.dtmax);
    outputEnergy[0].println("Length: "+fs.Length+" Mass: "+fs.Mass+" Points: "+fs.numOfpoints+" Stiffness: "+fs.stiffness+" Damping: "+fs.damping+" dt: "+fs.dtmax);
  }
  void InitSheetOut( FlexibleSheet [] fs ) {
    for ( int i=0; i<fs.length; i++ ) {
      outputPos[i].println("Length: "+fs[i].Length+" Mass: "+fs[i].Mass+" Points: "+fs[i].numOfpoints+" Stiffness: "+fs[i].stiffness+" Damping: "+fs[i].damping+" dt: "+fs[i].dtmax);
      outputVel[i].println("Length: "+fs[i].Length+" Mass: "+fs[i].Mass+" Points: "+fs[i].numOfpoints+" Stiffness: "+fs[i].stiffness+" Damping: "+fs[i].damping+" dt: "+fs[i].dtmax);
      outputForce[i].println("Length: "+fs[i].Length+" Mass: "+fs[i].Mass+" Points: "+fs[i].numOfpoints+" Stiffness: "+fs[i].stiffness+" Damping: "+fs[i].damping+" dt: "+fs[i].dtmax);
      outputEnergy[i].println("Length: "+fs[i].Length+" Mass: "+fs[i].Mass+" Points: "+fs[i].numOfpoints+" Stiffness: "+fs[i].stiffness+" Damping: "+fs[i].damping+" dt: "+fs[i].dtmax);
    }
  }

  void InitPoints( PrintWriter out ) {
    out.println("Points: "+Ncp);
  }

  void saveInfoCPoints( float t, int d ) {
    outputPos[d].println("============= t = "+t+" ================");
    outputVel[d].println("============= t = "+t+" ================");
    outputForce[d].println("============= t = "+t+" ================");
    outputEnergy[d].println("============= t = "+t+" ================");

    for (int i=0; i<myCPoints[d].size(); i++) {
      ControlPoint cp = myCPoints[d].get(i);
      outputPos[d].println(cp.position.x + " " + cp.position.y);
      outputVel[d].println(cp.velocity.x + " " + cp.velocity.y);
      outputForce[d].println(cp.force.x + " " + cp.force.y);
      float EE = 0.5*cp.mass*sq(cp.velocity.mag());
      outputEnergy[d].println(EE);
    }
  }
  void saveInfoCPoints( float t ) { saveInfoCPoints( t, 0 ); }
  void saveInfoCPoints() { saveInfoCPoints( 0, 0 ); }
  
  void saveInfoCPointsSheet( float t, int d ) {
    outputPos[d].println("============= t = "+t+" ================");
    outputVel[d].println("============= t = "+t+" ================");
    outputForce[d].println("============= t = "+t+" ================");

    for (int i=0; i<myCPoints[d].size(); i++) {
      ControlPoint cp = myCPoints[d].get(i);
      outputPos[d].println(cp.position.x + " " + cp.position.y);
      outputVel[d].println(cp.velocity.x + " " + cp.velocity.y);
      outputForce[d].println(cp.force.x + " " + cp.force.y);
    }
  }
  void saveInfoCPointsSheet( float t ) { saveInfoCPoints( t, 0 ); }
  void saveInfoCPointsSheet() { saveInfoCPoints( 0, 0 ); }

  void saveInfoSheet( float t, PVector g, float b, int d ) {
    float gMag = g.mag();
    saveInfoCPointsSheet( t, d );
    float EE = 0;

    outputEnergy[d].println("============= t = "+t+" ================");
    for (int i=0; i<mySprings[d].size(); i++) {
      Spring spr = mySprings[d].get(i);
      EE += .5 * spr.stiffness * spr.getStretch() * spr.getStretch();
    }
    for (int i=0; i<myCPoints[d].size(); i++) {
      ControlPoint cpoi = myCPoints[d].get(i);
      EE += .5 * cpoi.mass * cpoi.velocity.mag() * cpoi.velocity.mag();
      EE += cpoi.mass * gMag * (b - cpoi.position.y);
    }
    outputEnergy[d].println(EE);
  }
  void saveInfoSheet( float t, int d ) { saveInfoSheet( t, new PVector(0,0), 0, d); }
  void saveInfoSheet( float t ) { saveInfoSheet( t, new PVector(0,0), 0, 0); }
  void saveInfoSheet() { saveInfoSheet( 0, new PVector(0,0), 0, 0); }

  void closeInfos() {
    if (posFlag) terminateFile(outputPos);
    if (velFlag) terminateFile(outputVel);
    if (forFlag) terminateFile(outputForce);
    if (eneFlag) terminateFile(outputEnergy);
  }

  void terminateFile(PrintWriter [] file) {
    for (int j=0; j<file.length; j++) {
      file[j].flush(); // Writes the remaining data to the file
      file[j].close(); // Finishes the file
    }
  }

} // end of class