/**********************************************************************
      WriteInfo class: Produces the output txt files containing 
      information on the test, the sheet(s) and the control 
      points. Total energy and length of the sheet are recorded,
      as well as position, velocity and force on each of the 
      control points.

**********************************************************************/
class WriteInfo {
  
  int Nfs; 
  FlexibleSheet [] sheets;
  PrintWriter [][] cpointsFiles;
  PrintWriter [] energyFiles;
  PrintWriter [] genInfoFiles;

  //================= Constructor ====================//
  WriteInfo( FlexibleSheet fs ) {
    Nfs = 1;
    int Ncp = fs.numOfpoints;
    sheets = new FlexibleSheet[Nfs];
    sheets[0] = fs;
    
    cpointsFiles = new PrintWriter[Nfs][Ncp];
    for (int i=0; i<Ncp; i++) {
      cpointsFiles[0][i] = createWriter("./info/sheet0/cpoints"+i+".txt");
    }
    energyFiles = new PrintWriter[Nfs];
    energyFiles[0] = createWriter("./info/sheet0/energy.txt");
    genInfoFiles = new PrintWriter[Nfs];
    genInfoFiles[0] = createWriter("./info/sheet0/generalInfo.txt");
    this.writeGenInfo();
  }

  WriteInfo(FlexibleSheet [] fs) {
    Nfs = fs.length;
    sheets = new FlexibleSheet[Nfs];
    for (int i=0; i<Nfs; i++) sheets[i] = fs[i];
    
    int Ncp = fs[0].numOfpoints;
    for (FlexibleSheet ff : fs) {
      if (ff.numOfpoints>Ncp) Ncp = ff.numOfpoints;
    }
    
    cpointsFiles = new PrintWriter[Nfs][Ncp];
    for (int j=0; j<Nfs; j++) {
      for (int i=0; i<fs[j].numOfpoints; i++) {
        cpointsFiles[j][i] = createWriter("./info/sheet"+j+"/cpoints"+i+".txt"); 
      }
    }
    energyFiles = new PrintWriter[Nfs];
    for (int j=0; j<Nfs; j++) energyFiles[j] = createWriter("./info/sheet"+j+"/energy.txt");
    genInfoFiles = new PrintWriter[Nfs];
    for (int j=0; j<Nfs; j++) genInfoFiles[j] = createWriter("./info/sheet"+j+"/generalInfo.txt");
    this.writeGenInfo();
  }

  //================= Methods ====================//
  // Write some general information on the sheet
  void writeGenInfo() {
    int N = genInfoFiles.length;
    for (int i=0; i<N; i++) {
      genInfoFiles[i].println("Length: "+sheets[i].Length);
      genInfoFiles[i].println("Mass: "+sheets[i].Mass);
      genInfoFiles[i].println("Points: "+sheets[i].numOfpoints);
      genInfoFiles[i].println("Stiffness: "+sheets[i].stiffness);
      genInfoFiles[i].println("Damping: "+sheets[i].damping);
      genInfoFiles[i].println("dt: "+sheets[i].dtmax);
      genInfoFiles[i].println("----------------------");
    }
  }
  
  // Add extra lines in the general info file (should be object-specific)
  void addGenInfo( int id, String str ) {
    genInfoFiles[id].println( str );
  }
  
  // Calculate and write the total energy and the length of the sheet
  void writeEnergy( float t, boolean fl ) {
    for (int i=0; i<Nfs; i++) {
      FlexibleSheet fs = sheets[i];
      float EE = 0;
      for (int j=0; j<fs.numOfsprings; j++) {
        Spring spr = fs.springs[j];
        EE += .5 * spr.stiffness * spr.getStretch() * spr.getStretch();
      }
      for (int j=0; j<fs.numOfpoints; j++) {
        ControlPoint cpoi = fs.cpoints[j];
        EE += .5 * cpoi.mass * cpoi.velocity.mag() * cpoi.velocity.mag();
        if (fl) EE += cpoi.mass * (fs.window.idy(fs.window.dy) - cpoi.position.y);
      }
    energyFiles[i].println(t+","+EE+","+sheets[i].CurrentLength());
    }
  }
  
  // Write all info of control points 
  void writeCPoints( float t ) {
    for (int i=0; i<Nfs; i++) {
      for (int j=0; j<sheets[i].numOfpoints; j++) {
        sheets[i].cpoints[j].dampInfo( cpointsFiles[i][j], t );
      }
    }
  }
  
  // Method to call externally, handles all output
  void dampAllInfo( float t, boolean fl ) {
    this.writeEnergy( t, fl );
    this.writeCPoints( t );
  }
  
  // Gracefully terminate the writing of the files.
  void terminateFiles() {
    for (int i=0; i<Nfs; i++) {
      energyFiles[i].flush();
      energyFiles[i].close();
      genInfoFiles[i].flush();
      genInfoFiles[i].close();
      for (int j=0; j<sheets[i].numOfpoints; j++) {
        cpointsFiles[i][j].flush();
        cpointsFiles[i][j].close();
      }
    }
  }

} // end of class
