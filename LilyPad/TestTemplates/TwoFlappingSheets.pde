Window view; // window in h-units
FloodPlot plot;

boolean saveimg = false;

int nx = (int)pow(2,6); // x-dir
int ny = (int)pow(2,6); // y-dir

float Re = 100;
float L = nx/4.;
float thick = 1;
int resol = 1;
float M = 5;
float stiff = 100;
PVector lpos = new PVector( nx/3., 0.5*ny );
PVector align = new PVector(1, 0);
float separ = L/2.;

float t=0;
float dt;

BDIM flow;
FlexibleSheet [] sheet;
BodyUnion bodies;
CollisionHandler collide;
WriteInfo writer;

void settings(){
    size(600, 600);
}

void setup() {
  view = new Window( 1, 1, nx, ny, 0, 0, width, height);
  sheet = new FlexibleSheet[2];
  
  sheet[0] = new FlexibleSheet(L, thick, M, resol, stiff, new PVector(lpos.x, lpos.y-separ/2.), align, view);
  sheet[0].cpoints[0].makeFixed(); // pinning leading point
  sheet[1] = new FlexibleSheet(L, thick, M, resol, stiff, new PVector(lpos.x, lpos.y+separ/2.), align, view);
  sheet[1].cpoints[0].makeFixed(); // pinning leading point
  
  plot = new FloodPlot(view); // standard window
  
  dt = sheet[0].dtmax;
  bodies = new BodyUnion(sheet[0], sheet[1]);
  
  flow = new BDIM(nx, ny, dt, bodies, 1/Re, true);
  
  plot.range = new Scale(-1,1);
  plot.hue = new Scale(100, 40);
  plot.setLegend("pressure");
  
  collide = new CollisionHandler( sheet );
  writer = new WriteInfo( sheet );
} // end of setup


void draw() {
  
  for (FlexibleSheet fs : sheet) fs.updateAlt(dt, flow);
  flow.update(bodies);
  
  for (FlexibleSheet fs : sheet) fs.updateAlt2(dt, flow);
  flow.update2();
  
  collide.HandleCollisions();
  
  plot.display(flow.p);
  for (FlexibleSheet fs : sheet) fs.display();
  
  writer.dampAllInfo( t, false );
  if (saveimg) saveFrame("movie/frame_######.png");
  t += dt;
}

void keyPressed() {
  writer.terminateFiles();
}


