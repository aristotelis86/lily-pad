/*********************************************************
                  Main Window! 

Click the "Run" button to Run the simulation.

Change the geometry, flow conditions, numercial parameters
 visualizations and measurments from this window. 

*********************************************************/
BDIM flow;
EllipseBody body;
FloodPlot flood;

void setup(){
  int n=(int)pow(2,6)+2, z=4; // number of grid points and display zoom (pixels/point)
  size(z*n,z*n);              // display window size

  body = new EllipseBody(n/3,n/2,n/8,new Window(n,n)); // define geom (EllipseBody, initial location, radius, 
  flow = new BDIM(n,n,1.5,body);
  flood = new FloodPlot(new Window(n,n));
  flood.range = new Scale(-.75,.75);
  flood.setLegend("vorticity");
}
void draw(){
  body.update();
  flow.update(body);
  flow.update2(body);
  flood.display(flow.u.vorticity());
  body.display();
}
void mousePressed(){body.mousePressed();}
void mouseReleased(){body.mouseReleased();}

