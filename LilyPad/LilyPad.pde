/*********************************************************
                  Main Window!

Click the "Run" button to Run the simulation.

Change the geometry, flow conditions, numercial parameters
visualizations and measurments from this window.

This screen has an example. Other examples are found at 
the top of each tab. Copy/paste them here to run, but you 
can only have one setup & run at a time.

*********************************************************/

Window view, view2; // window in h-units

boolean saveimg = true;

int nx = (int)pow(2,7); // x-dir
int ny = (int)pow(2,7); // y-dir

float L = nx/6.;
float M = 100;
float stiff = 20;

float xpos = nx/4.;
float ypos = ny/2.;
PVector align = new PVector(1, 0);

float t=0;
float dt;

float sinAmp = L/5.;
float sinN = 1.;

BDIM flow, flow2;
FlexibleSheet sheet,sheet2;
FloodPlot plot, plot2;


void setup() {
  size(800, 600);
  //view = new Window(nx, ny);
  view = new Window( 1, 1, nx, ny, 0, 0, 800, 300);
  view2 = new Window( 1, 1, nx, ny, 0, 300, 800, 300);
  
  
  sheet = new FlexibleSheet(L, M, stiff, xpos, ypos, align, view);
  sheet2 = new FlexibleSheet(L, M, stiff, xpos, ypos, align, view2);
  
  plot = new FloodPlot(view); // standard window
  plot2 = new FloodPlot(view2); // standard window
  
  sheet.cpoints[0].makeFixed();
  sheet2.cpoints[0].makeFixed();
  
  dt = sheet.dtmax;
  
  flow = new BDIM(nx, ny, dt, sheet, 0.01, true);
  flow2 = new BDIM(nx, ny, dt, sheet2, 0.01, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
  plot2.range = new Scale(-1,1);
  plot2.hue = new Scale(5, 220);
  plot2.setLegend("vorticity");
  
} // end of setup


void draw() {
  
  sheet.update(dt, flow);
  sheet.update2(dt, flow);
  
  flow.update(sheet);
  flow.update2();
  
  //flow.u.curl().display(-0.5,0.5);
  plot.display(flow.p);
  
  sheet.display();
  // -------------------- // 
  sheet2.update(dt, flow2);
  sheet2.update2(dt, flow2);
  
  flow2.update(sheet2);
  flow2.update2();
  
  //flow.u.curl().display(-0.5,0.5);
  plot2.display(flow2.u.curl());
  
  sheet2.display();

  t += dt;
  
  if (saveimg) saveFrame("movie/frame_######.png");
  
  //noLoop();
}
void mousePressed(){sheet.mousePressed();}    // user mouse...
void mouseReleased(){sheet.mouseReleased();}  // interaction methods
void mouseWheel(MouseEvent event){sheet.mouseWheel(event);}



void updateSheet(float t) {
  
  int nn = sheet.cpoints.length;
  float [] x = new float[nn];
  float [] y = new float[nn];
  float [] vx = new float[nn];
  float [] vy = new float[nn];
  
  x[0] = xpos;
  y[0] = ypos;
  vx[0] = 0.;
  vy[0] = 0.;
  for (int i = 1; i < nn; i++) {
    x[i] = (L/(nn-1)) + x[i-1];
    y[i] = (sinAmp * sin(sinN*PI*(x[i]-xpos)/L))*sin(2*t) + ypos;
    vx[i] = 0.;
    vy[i] = 2*cos(2*t)*(sinAmp * sin(sinN*PI*(x[i]-xpos)/L));
  }
  sheet.UpdateState(x, y, vx, vy);
  
}



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