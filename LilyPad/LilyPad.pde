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
int nx = (int)pow(2,8); // x-dir
int ny = (int)pow(2,8); // y-dir

float L = nx/4.;
float M = 10;
float stiff = 100;

float xpos = nx/4.;
float ypos = ny/2.;
PVector align = new PVector(1, 0);

float t=0;
float dt = 0.1;

float sinAmp = L/10.;
float sinN = 2.;

BDIM flow;
FlexibleSheet sheet;
//Body sheet;

void setup() {
  size(500, 500);
  view = new Window(nx, ny);
  
  sheet = new FlexibleSheet(L, M, stiff, xpos, ypos, align, view);
  
  // Create the distortion
  int nn = sheet.cpoints.length;
  float [] x = new float[nn];
  float [] y = new float[nn];
  
  x[0] = xpos;
  y[0] = ypos;
  for (int i = 1; i < nn; i++) {
    x[i] = (L/(nn-1)) + x[i-1];
    y[i] = sinAmp * sin(sinN*PI*(x[i]-xpos)/L) + ypos;
    sheet.cpoints[i].position.x = x[i];
    sheet.cpoints[i].position.y = y[i];
  }
  sheet.getOrth();
  dt = sheet.dtmax;
  
  flow = new BDIM(nx, ny, dt, sheet, 0.01, true);
  
} // end of setup


void draw() {
  
  //if (t<3) sheet.move();
  //sheet.follow();
  //sheet.update(dt, new PVector(0,0));
  flow.update(sheet);
  flow.update2();
  
  flow.u.curl().display(-0.5,0.5);
  sheet.display();

  t += dt;
}
void mousePressed(){sheet.mousePressed();}    // user mouse...
void mouseReleased(){sheet.mouseReleased();}  // interaction methods
void mouseWheel(MouseEvent event){sheet.mouseWheel(event);}





//// Circle that can be dragged by the mouse
//BDIM flow;
//Body body;
//FloodPlot flood;

//void setup(){
//  size(700,700);                             // display window size
//  int n=(int)pow(2,7);                       // number of grid points
//  float L = n/8.;                            // length-scale in grid units
//  Window view = new Window(n,n);

//  body = new CircleBody(n/3,n/2,L,view);     // define geom
//  flow = new BDIM(n,n,1.5,body);             // solve for flow using BDIM
//  flood = new FloodPlot(view);               // intialize a flood plot...
//  flood.setLegend("vorticity",-.5,.5);       //    and its legend
//}
//void draw(){
//  body.follow();                             // update the body
//  flow.update(body); flow.update2();         // 2-step fluid update
//  flood.display(flow.u.curl());              // compute and display vorticity
//  body.display();                            // display the body
//}
//void mousePressed(){body.mousePressed();}    // user mouse...
//void mouseReleased(){body.mouseReleased();}  // interaction methods
//void mouseWheel(MouseEvent event){body.mouseWheel(event);}