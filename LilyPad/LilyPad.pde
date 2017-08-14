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

float L = nx/6.;
float M = 100;
float stiff = 50;

float xpos = nx/4.;
float ypos = ny/2.;
PVector align = new PVector(1, 0);

float t=0;
float dt;

float sinAmp = L/5.;
float sinN = 1.;

BDIM flow;
FlexibleSheet sheet;
FloodPlot plot;
Field press;

void setup() {
  size(800, 800);
  view = new Window(nx, ny);
  
  sheet = new FlexibleSheet(L, M, stiff, xpos, ypos, align, view);
  plot = new FloodPlot(view); // standard window
  
  press = new Field(nx, ny, 0, 5);
  for( int i=0; i<nx; i++){
  for( int j=0; j<ny; j++){
    press.a[i][j] = 24.*(0.5-j/(float)(ny-1));
    //press.a[i][j] = 24.*(0.5-i/(float)(nx-1));
  }}
  
  // Create the distortion
  sheet.cpoints[0].makeFixed();
  //int nn = sheet.cpoints.length;
  //sheet.cpoints[nn-1].makeFixed();
  
  dt = sheet.dtmax;
  
  flow = new BDIM(nx, ny, dt, sheet, 0.01, true);
  
  //PVector [] sheetPf;
  
  //sheetPf = sheet.pressForcePoints ( press );
  //for (int j=0; j<sheetPf.length; j++) println(sheetPf[j]);
  //println(press.linear(sheet.cpoints[8].position.x, sheet.cpoints[8].position.y));
  
  plot.range = new Scale(-2,2);
  plot.hue = new Scale(200, 140);
  
  //println(press.a[38][65]);
  //PVector [] printForce;
  //printForce = sheet.pressForcePoints(press);
  //for (int i=0; i<nx/6; i++) println(printForce[i]);
  
} // end of setup


void draw() {
  
  //if (t<3) sheet.move();
  //sheet.follow();
  //updateSheet(t);
  
  //sheet.update(dt, press);
  //sheet.update2(dt, press);
  sheet.update(dt, flow.p);
  sheet.update2(dt, flow.p);
  
  flow.update(sheet);
  flow.update2();
  
  //flow.u.curl().display(-0.5,0.5);
  plot.display(flow.p);
  //plot.display(press);
  sheet.display();

  t += dt;
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