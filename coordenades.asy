import geometry;
import graph;

picture Quadricula(real x0,
             		   real x1,
             		   real y0,
             		   real y1,
             		   pen color=gray(0.8),
             		   int divs = 4,
             		   real unit = 1.0,
                   real scaled=1.0)
{
  pen e = color + linewidth(scaled*0.15mm);
  pen t = dotted + color + linewidth(scaled*0.06mm);
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
                    real tickx=0.15,
                    real ticky=0.15,
                    int divs=4,
                    int sizex=1,
                    int sizey=sizex,
                    bool base = false,
                    bool quadricula=true,
                    real scaled=1.0)
{
  picture pic;

  if(quadricula)
    add(pic,Quadricula(x0,x1,y0,y1,divs=divs,scaled=scaled));
  pen e = color + linewidth(scaled*0.25mm);
  pen t = color + linewidth(scaled*0.15mm);
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
      draw(pic,(i,0)--(i,-tickx),e);
      label(pic,scale(scaled)*format("$%i$",i),(i,-tickx*scaled),S,e);
      for(int k=1;k<divs;++k)
        draw(pic,(i+signe*step*k,0)--(i+signe*step*k,-tickx*scaled/2),t);
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
      draw(pic,(0,i)--(-ticky,i),e);
      label(pic,scale(scaled)*format("$%i$",i),(-ticky*scaled,i),W,e);
      for(int k=1;k<divs;++k)
         draw(pic,(0,i+signe*step*k)--(-ticky/2,i+signe*step*k),t);
    }
  }

  if (base)
  {
    pen p = color + linewidth(scaled*0.55mm);
    draw(pic,(0,0)--(1,0),p,Arrow(scaled*2mm));
    draw(pic,(0,0)--(0,1),p,Arrow(scaled*2mm));
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
                      bool quadricula=true,
                      real scaled=1.0)
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
    draw(pic,(i * pi/2,0)--(i * pi/2,-0.15*scaled),e);
    for(int k=1;k<divs;++k)
      draw((i*pi/2+signe*step*k,0)--(i*pi/2+signe*step*k,-0.075*scaled),t);

    if(i == 1)
    {
      label(pic,scale(scaled*1.09)*"$\frac{\pi}{2}$",(i * pi/2,-0.15*scaled),S);
      continue;
    }
    if(i == -1)
    {
      label(pic,scale(scaled*1.09)*"$-\frac{\pi}{2}$",(i * pi/2,-0.15*scaled),S);
      continue;
    }
    if (i == 2)
    {
      label(pic,scale(scaled)*"$\pi$",(i * pi/2,-0.15*scaled),S);
      continue;
    }
    if (i == -2)
    {
      label(pic,scale(scaled)*"$-\pi$",(i * pi/2,-0.15*scaled),S);
      continue;
    }
    if (i % 2 == 0)
    {
      label(pic,scale(scaled)*format("$%i\pi$",i # 2),(i * pi/2,-0.15*scaled),S);
      continue;
    }
    if (i > 0)
      label(pic,scale(scaled*1.09)*format("$\frac{%i\pi}{2}$",i),(i * pi/2,-0.15*scaled),S);
    else
      label(pic,scale(scaled*1.09)*format("$-\frac{%i\pi}{2}$",-i),(i * pi/2,-0.15*scaled),S);
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
    draw(pic,(0,i*sizey)--(-0.15*scaled,i*sizey),e);
    label(pic,scale(scaled)*format("$%i$",i),(-0.15*scaled,i*sizey),W,e);
    for(int k=1;k<divs;++k)
      draw(pic,(0,i*sizey+signe*step*k)--(-0.075*scaled,i*sizey+signe*step*k),t);
  }

  return pic;
}

void Canonica(int x0,
              int x1,
              int y0,
              int y1,
              real scaled=1.0)
{
  picture e = Coordenades(x0,x1,y0,y1,scaled=scaled);
  add(scale(scaled) * e);
}

transform Referencia(pair origen=(0,0),
                     pair u1=(1,0),
                     pair u2=(0,1)
                    )
{
  transform tr = (origen.x,origen.y,u1.x,u2.x,u1.y,u2.y);
  return tr;
}

transform Base(pair u1=(1,0),
               pair u2=(0,1)
              )
{
  transform tr = (0,0,u1.x,u2.x,u1.y,u2.y);
  return tr;
}

transform TransformacioLineal(real[][] A)
{
  transform tr = (0,0,A[0][0],A[0][1],A[1][0],A[1][1]);
  return tr;
}

transform TransformacioAfi(real[] B,real[][] A)
{
  transform tr = (B[0],B[1],A[0][0],A[0][1],A[1][0],A[1][1]);
  return tr;
}

void ElipseSimple(pair o,
                  real a2,
                  real b2,
                  int x=10,
                  int y=8,
                  real scaled=1.0)
{
  picture rp = Coordenades(-x,x,-y,y,color=red,base=true,quadricula=false,scaled=scaled);
  real a, c;
  pair F;
  if(a2 > b2)
  {
    a = sqrt(a2);
    c = sqrt(a2 - b2);
    F = (c,0);
  }
  else
  {
    a = sqrt(b2);
    c = sqrt(b2 - a2);
    F = (0,c);
  }
  ellipse el =  ellipse(F,-F,a);
  draw(rp,el,blue+linewidth(scaled*0.55mm));
	dot(rp,F,green+scaled*2mm);
	dot(rp,-F,green+scaled*2mm);
  transform tr = (o.x,o.y,1,0,0,1);
  path cl = (-x,-y)--(x,-y)--(x,y)--(-x,y)--cycle;
  clip(rp,cl);
  add(scale(scaled) * tr * rp);
}

void Elipse(pair o,
            pair u,
            real a2,
            real b2,
            int x=10,
            int y=8,
            real scaled=1.0)
{
  if(dot(u,(0,1)) == 0)
  {
  	ElipseSimple(o,a2,b2,x,y,scaled=scaled);
  	return;
  }
  if(dot(u,(1,0)) == 0)
  {
  	ElipseSimple(o,b2,a2,y,x,scaled=scaled);
  	return;
  }
  picture rp = Coordenades(-x,x,-y,y,color=red,base=true,quadricula=false,scaled=scaled);
  real a, c;
  pair F;
  if(a2 > b2)
  {
    a = sqrt(a2);
    c = sqrt(a2 - b2);
    F = (c,0);
  }
  else
  {
    a = sqrt(b2);
    c = sqrt(b2 - a2);
    F = (0,c);
  }
  ellipse el =  ellipse(F,-F,a);
  draw(rp,el,blue+linewidth(scaled*0.55mm));
	dot(rp,F,green+scaled*2mm);
	dot(rp,-F,green+scaled*2mm);
  real l = length(u);
  transform tr = (o.x,o.y,u.x/l,-u.y/l,u.y/l,u.x/l);
  path cl = (-x,-y)--(x,-y)--(x,y)--(-x,y)--cycle;
  clip(rp,cl);
  add(tr * rp);
}

void HiperbolaSimple(pair o,
                     real a2,
                     real b2,
                     int x=10,
                     int y=8,
                     bool xdir,
                     real scaled=1.0)
{
  pen lx;
  picture rp;
  lx = green+linewidth(scaled*0.4mm);
  real a = sqrt(a2);
  real b = sqrt(b2);
  transform tr = (o.x,o.y,1,0,0,1);
  transform th;
  if(xdir)
  {
    rp = Coordenades(-x,x,-y,y,color=red,base=true,quadricula=false,scaled=scaled);
    th = (0,0,1,0,0,1);
  }
  else
  {
    rp = Coordenades(-y,y,-x,x,color=red,base=true,quadricula=false,scaled=scaled);
    th = (0,0,0,1,-1,0);
  }
  hyperbola h = hyperbola((0,0),a,b);
  real m = b/a;
  real c = sqrt(a^2+b^2);
  picture rh;
  draw(rh,(-x,-x*m)--(x,x*m),lx);
  draw(rh,(-x,x*m)--(x,-x*m),lx);
  dot(rh,(c,0),green+scaled*2mm);
  dot(rh,(-c,0),green+scaled*2mm);
  draw(rh,h,blue+linewidth(scaled*0.55mm));
  real l = 1.0*x;
  if (y > x)
    l = 1.0*y;
  path cl = circle((0,0),l);
  clip(rp,cl);
  clip(rh,cl);
  add(scale(scaled) * tr * rp);
  add(scale(scaled) * tr * th * rh);
}

void Hiperbola(pair o,
               pair u,
               real a2,
               real b2,
               int x=10,
               int y=8,
               real scaled=1.0)
{
  if(dot(u,(0,1)) == 0)
  {
  	HiperbolaSimple(o,a2,b2,xdir=true,scaled=scaled);
  	return;
  }
  if(dot(u,(1,0)) == 0)
  {
  	HiperbolaSimple(o,a2,b2,xdir=false,scaled=scaled);
  	return;
  }
  pen lx;
  lx = green+linewidth(scaled*0.4mm);
  picture rp = Coordenades(-x,x,-y,y,color=red,base=true,quadricula=false,scaled=scaled);
	real a = sqrt(a2);
	real b = sqrt(b2);
  hyperbola h = hyperbola((0,0),a,b);
  real m = b/a;
	real c = sqrt(a^2+b^2);
  draw(rp,(-x,-x*m)--(x,x*m),lx);
  draw(rp,(-x,x*m)--(x,-x*m),lx);
	dot(rp,(c,0),green+scaled*2mm);
	dot(rp,(-c,0),green+scaled*2mm);
  draw(rp,h,blue+linewidth(scaled*0.55mm));
  real l = length(u);
  transform tr = (o.x,o.y,u.x/l,-u.y/l,u.y/l,u.x/l);
	l = x;
	if (y > x)
		l = y;
  path cl = circle((0,0),l);
  clip(rp,cl);
  add(tr * rp);
}

void ParabolaSimple(pair v,
              	   pair f,
              	   int x=10,
                   int y=8,
                   real scaled=1.0)
{
  pen lx;
  pen lx = green+linewidth(scaled*0.4mm);
  picture rp = Coordenades(-x,x,-y,y,color=red,base=true,quadricula=false,scaled=scaled);
  pair s = f - v;
  real par = 2 * length(s);
  transform t = (v.x,v.y,1,0,0,1);
  parabola pa = parabola(inverse(t)*f,(0,0));
  draw(rp,pa,blue+linewidth(scaled*0.55mm));
  dot(rp,inverse(t)*f,green+scaled*2mm);
  real dp1 = dot(s,(1,0));
  real dp2 = dot(s,(0,1));
  if(dp1 == 0)
  {
    if(dp2 < 0)
      par = -par;
  	draw(rp,(-x,-par/2)--(x,-par/2),lx);
  }
  else
  {
    if(dp1 < 0)
      par = -par;
    draw(rp,(-par/2,-y)--(-par/2,y),lx);
  }
  path cl = (-x,-y)--(x,-y)--(x,y)--(-x,y)--cycle;
  clip(rp,cl);
  add(t * rp);
}

void Parabola(pair v,
              pair f,
              int d=1,
              int x=10,
              int y=8,
              real scaled=1.0)
{
  pair s = f - v;
  if(dot(s,(1,0)) == 0 || dot(s,(0,1)) == 0)
  {
  	ParabolaSimple(v,f,x,y,scaled=scaled);
  	return;
  }
  pen lx = green+linewidth(scaled*0.4mm);
  picture rp = Coordenades(-x,x,-y,y,color=red,base=true,quadricula=false,scaled=scaled);
  pair s = f - v;
  real par = 2 * length(s);
  if(d < 0)
  {
  	s = -s;
  	par = -par;
  }
  pair u = (s.y,-s.x);
  parabola pa = parabola((0,par/2),(0,0));
  draw(rp,(-x,-par/2)--(x,-par/2),lx);
  dot(rp,(0,par/2),green+scaled*2mm);
  draw(rp,pa,blue+linewidth(scaled*0.55mm));
  real l = length(u);
  transform tr = (v.x,v.y,u.x/l,-u.y/l,u.y/l,u.x/l);
  path cl = (-x,-y)--(x,-y)--(x,y)--(-x,y)--cycle;
  clip(rp,cl);
  add(tr * rp);
}
