import geometry;
import graph;

picture Quadricula(real x0,
             	   real x1,
             	   real y0,
             	   real y1,
             	   pen color=gray(0.8),
             	   int divs = 4,
             	   real unit = 1.0)
{
  pen e = color + linewidth(0.15mm);
  pen t = dotted + color + linewidth(0.06mm);
  real step = unit / divs;
  picture pic;
  for(real i=x0;i<=x1;i+=unit)
  {
    draw(pic,(i,y0)--(i,y1),e);
    if(i < x1)
      for(int k=1;k<divs;++k)
        draw(pic,(i+step*k,y0)--(i+step*k,y1),t);

  }
  for(real i=y0;i<=y1;i+=unit)
  {
    draw(pic,(x0,i)--(x1,i),e);
    if(i < y1)
      for(int k=1;k<=divs;++k)
        if (i+0.25*k*unit <= y1)
         draw(pic,(x0,i+step*k)--(x1,i+step*k),t);
  }
  return pic;
}

picture Coordenades(int x0,
                    int x1,
                    int y0,
                    int y1,
                    pen color=black,
		    real scaletext=0.5,
                    int divs=4,
                    int sizex=1,
                    int sizey=sizex,
                    bool base=false,
                    bool quadricula=true)
{
  picture pic;

  if(quadricula)
    add(pic,Quadricula(x0,x1,y0,y1,divs=divs));
  pen e = color + linewidth(0.25mm);
  pen t = color + linewidth(0.15mm);
  draw(pic,(x0,0)--(x1,0),e);
  draw(pic,(0,y0)--(0,y1),e);

  for(int i=x0+1;i<=x1;i+=sizex)
  {
    real step = sizex / divs;
    int signe = 1;
    if(i > 0)
      signe = -1;
    if (i == 0)
      continue;
    if(i < x1)
    {
      draw(pic,(i,0)--(i,-0.15),e);
      label(pic,scale(scaletext)*format("$%i$",i),(i,-0.15),S,e);
      for(int k=1;k<divs;++k)
        draw(pic,(i+signe*step*k,0)--(i+signe*step*k,-0.075),t);
    }
  }

  for(int i=y0+1;i<=y1;i+=sizey)
  {
    real step = sizey / divs;
    int signe = 1;
    if(i > 0)
      signe = -1;
    if(i == 0)
      continue;
    if(i < y1)
    {
      draw(pic,(0,i)--(-0.15,i),e);
      label(pic,scale(scaletext)*format("$%i$",i),(-0.15,i),W,e);
      for(int k=1;k<divs;++k)
         draw(pic,(0,i+signe*step*k)--(-0.075,i+signe*step*k),t);
    }
  }

  if (base)
  {
    pen p = color + linewidth(0.35mm);
    draw(pic,(0,0)--(1,0),p,Arrow(1.5mm));
    draw(pic,(0,0)--(0,1),p,Arrow(1.5mm));
  }
  return pic;
}

picture PiCoordenades(int p0,
                      int p1,
                      int y0,
                      int y1,
                      pen color=black,
                      int divs=4,
                      int sizex=1,
                      int sizey=sizex,
                      bool quadricula=true)
{
  picture pic;
  real x0 = p0 * pi,
       x1 = p1 * pi;

  if(quadricula)
    add(pic,Quadricula(floor(x0),ceil(x1),sizey*y0-1,sizey*y1+1,divs=divs));
  pen e = color + linewidth(0.25mm);
  pen t = color + linewidth(0.15mm);
  draw(pic,(floor(x0),0)--(ceil(x1),0),e);
  draw(pic,(0,y0*sizey-1)--(0,y1*sizey+1),e);

  for(int i=2p0;i <= 2p1;++i)
  {
    int signe = 1;
    if(i > 0)
      signe = -1;
    real step = pi / (2*divs);
    if(i == 0)
      continue;
    draw(pic,(i * pi/2,0)--(i * pi/2,-0.15),e);
    for(int k=1;k<divs;++k)
      draw((i*pi/2+signe*step*k,0)--(i*pi/2+signe*step*k,-0.075),t);

    if(i == 1)
    {
      label(pic,scale(1.09)*"$\frac{\pi}{2}$",(i * pi/2,-0.15),S);
      continue;
    }
    if(i == -1)
    {
      label(pic,scale(1.09)*"$-\frac{\pi}{2}$",(i * pi/2,-0.15),S);
      continue;
    }
    if (i == 2)
    {
      label(pic,"$\pi$",(i * pi/2,-0.15),S);
      continue;
    }
    if (i == -2)
    {
      label(pic,"$-\pi$",(i * pi/2,-0.15),S);
      continue;
    }
    if (i % 2 == 0)
    {
      label(pic,format("$%i\pi$",i # 2),(i * pi/2,-0.15),S);
      continue;
    }
    if (i > 0)
      label(pic,scale(1.09)*format("$\frac{%i\pi}{2}$",i),(i * pi/2,-0.15),S);
    else
      label(pic,scale(1.09)*format("$-\frac{%i\pi}{2}$",-i),(i * pi/2,-0.15),S);
  }

  divs *= sizey;
  for(int i=y0;i<=y1;i+=1)
  {
    real step = sizey / divs;
    int signe = 1;
    if(i > 0)
      signe = -1;
    if(i == 0)
      continue;
    draw(pic,(0,i*sizey)--(-0.15,i*sizey),e);
    label(pic,format("$%i$",i),(-0.15,i*sizey),W,e);
    for(int k=1;k<divs;++k)
      draw(pic,(0,i*sizey+signe*step*k)--(-0.075,i*sizey+signe*step*k),t);
  }

  return pic;
}

void Canonica(int x0,
              int x1,
              int y0,
              int y1,
	      real scaletext=0.5)
{
  picture e = Coordenades(x0,x1,y0,y1,scaletext=scaletext);
  add(e);
}

void Ellipse(pair o,pair u,real a2,real b2,int x=10,int y=8)
{
  picture rp = Coordenades(-x,x,-y,y,color=red,base=true,quadricula=false);
  if (b2 > a2)
    return;
  real a = sqrt(a2);
  real c = sqrt(a2 - b2);
  pair F = (c,0);
  ellipse el =  ellipse(F,-F,a);
  draw(rp,el,blue+linewidth(0.35mm));
  dot(rp,F,green+0.75mm);
  dot(rp,-F,green+0.75mm);
  real l = length(u);
  transform tr = (o.x,o.y,u.x/l,-u.y/l,u.y/l,u.x/l);
  path cl = (-x,-y)--(x,-y)--(x,y)--(-x,y)--cycle;
  clip(rp,cl);
  add(tr * rp);
}

void Hiperbola(pair o,pair u,real a2,real b2,int x=10,int y=8)
{
  pen lx;
  lx = green+linewidth(0.2mm);
  picture rp = Coordenades(-x,x,-y,y,color=red,base=true,quadricula=false);
  real a = sqrt(a2);
  real b = sqrt(b2);
  hyperbola h = hyperbola((0,0),a,b);
  real m = b/a;
	real c = sqrt(a^2+b^2);
  draw(rp,(-x,-x*m)--(x,x*m),lx);
  draw(rp,(-x,x*m)--(x,-x*m),lx);
  dot(rp,(c,0),green+0.75mm);
  dot(rp,(-c,0),green+0.75mm);
  draw(rp,h,blue+linewidth(0.35mm));
  real l = length(u);
  transform tr = (o.x,o.y,u.x/l,-u.y/l,u.y/l,u.x/l);
  l = x;
  if (y > x)
    l = y;
  path cl = circle((0,0),l);
  clip(rp,cl);
  add(tr * rp);
}

void Parabola(pair v,pair f,int x=10,int y=8)
{
  pen lx;
  pen lx = green+linewidth(0.2mm);
  picture rp = Coordenades(-x,x,-y,y,color=red,base=true,quadricula=false);
  pair s = f - v;
  pair u = (s.y,-s.x);
  real par = 2 * length(s);
  parabola pa = parabola((0,par/2),(0,0));
  draw(rp,(-x,-par/2)--(x,-par/2),lx);
  dot(rp,(0,par/2),green+0.75mm);
  draw(rp,pa,blue+linewidth(0.35mm));
  real l = length(u);
  transform tr = (v.x,v.y,u.x/l,-u.y/l,u.y/l,u.x/l);
  path cl = (-x,-y)--(x,-y)--(x,y)--(-x,y)--cycle;
  clip(rp,cl);
  add(tr * rp);
}
