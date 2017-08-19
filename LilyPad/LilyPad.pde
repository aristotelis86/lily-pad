/*********************************************************
                  Main Window!

Click the "Run" button to Run the simulation.

Change the geometry, flow conditions, numercial parameters
visualizations and measurments from this window.

This screen has an example. Other examples are found at 
the top of each tab. Copy/paste them here to run, but you 
can only have one setup & run at a time.

*********************************************************/

Window view; // window in h-units

int nx = (int)pow(2,7); // x-dir
int ny = (int)pow(2,7); // y-dir

VectorField testField;
FloodPlot plot;
Field u, v;

float L = nx/6.;
float M = 100;
float stiff = 20;

float xpos = nx/6.;
float ypos = ny/2.;
PVector align = new PVector(1, 0);

float t=0;
float dt;

FlexibleSheet sheet;

void setup() {
  size(800, 600);
  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  Field u = new Field(nx,ny);
  Field v = new Field(nx,ny);
  
  for( int i=0; i<nx; i++){
  for( int j=0; j<ny; j++){
    
    float myX = (i-nx/2.)/(float(nx)-1)-.1;
    float myY = (j-ny/2.)/(float(ny)-1);
    float myR = myX*myX + myY*myY;
    
    u.a[i][j] = (.1/myR)*cos(myY/myX);
    v.a[i][j] = (.1/myR)*sin(myY/myX);
    //u.a[i][j] = 2.5*j/(ny-1);
    //v.a[i][j] = 1.5*i/(nx-1);
  }}
  
  testField = new VectorField( u, v );
  
  plot = new FloodPlot(view);
  plot.range = new Scale(-5,5);
  plot.hue = new Scale(5, 220);
  plot.setLegend("vorticity");
  
  sheet = new FlexibleSheet(L, M, stiff, xpos, ypos, align, view);
  sheet.cpoints[0].makeFixed(); // pinning leading point
  dt = sheet.dtmax;
  
  
} // end of setup


void draw() {
  
  sheet.update(dt, testField, 1000);
  sheet.update2(dt, testField, 1000);
  
  plot.display(testField.curl());
  testField.display(1, 10);
  sheet.display();
  
  t += dt;
  //noLoop();
}




//Window view; // window in h-units

//boolean saveimg = false;

//int nx = (int)pow(2,7); // x-dir
//int ny = (int)pow(2,7); // y-dir

//float L = nx/6.;
//float M = 100;
//float stiff = 20;

//float xpos = nx/4.;
//float ypos = ny/2.;
//PVector align = new PVector(1, 0);

//float t=0;
//float dt;


//BDIM flow;
//FlexibleSheet sheet;
//FloodPlot plot;


//void setup() {
//  size(800, 600);
//  //view = new Window(nx, ny);
//  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  
//  sheet = new FlexibleSheet(L, M, stiff, xpos, ypos, align, view);
  
//  plot = new FloodPlot(view); // standard window
  
//  sheet.cpoints[0].makeFixed(); // pinning leading point
  
//  dt = sheet.dtmax;
  
//  flow = new BDIM(nx, ny, dt, sheet, 0.01, true);
  
//  plot.range = new Scale(-1,1);
//  plot.hue = new Scale(100, 40);
//  plot.setLegend("pressure");
  
//  //plot.range = new Scale(-1,1);
//  //plot.hue = new Scale(5, 220);
//  //plot.setLegend("vorticity");
  
//} // end of setup


//void draw() {
  
//  sheet.update(dt, flow);
//  sheet.update2(dt, flow);
  
//  flow.update(sheet);
//  flow.update2();
  
//  plot.display(flow.p);
//  //plot.display(flow2.u.curl());
  
//  sheet.display();
  
//  t += dt;
  
//  if (saveimg) saveFrame("movie/frame_######.png");
  
//  //noLoop();
//}