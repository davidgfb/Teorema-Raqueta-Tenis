// The MIT License
// Copyright © 2018 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// Cylinder intersection: https://www.shadertoy.com/view/4lcSRn
// Cylinder bounding box: https://www.shadertoy.com/view/MtcXRf
// Cylinder distance:     https://www.shadertoy.com/view/wdXGDr
// List of other 3D SDFs: https://www.shadertoy.com/playlist/43cXRl
// and https://iquilezles.org/articles/distfunctions
float sdCylinder(vec3 p, vec3 a, vec3 b, float r) {
    vec3  ba = b - a, pa = p - a;
    float baba = dot(ba,ba), paba = dot(pa,ba), 
        x = length(pa*baba-ba*paba) - r*baba, 
        y = abs(paba-baba/2.0)-baba/2.0, x2 = x*x, y2 = y*y*baba,
        d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    
    return sign(d)*sqrt(abs(d))/baba;
}

// vertical cylinder
float sdCylinder( vec3 p, float h, float r ) {
  vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(r,h);
  
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float map( in vec3 pos ) {
    float r = 0.4;
    vec3 x = vec3(1,0,0), y = vec3(0,1,0);
    
    return min(sdCylinder(pos, vec3(0), 2.0*y, r ),
        sdCylinder(pos, -x, x, r ));
}

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal( in vec3 pos ) {
    vec2 e = vec2(1.0,-1.0)*0.58;
    const float eps = 5e-4 ;
    return normalize( e.xyy*map( pos + e.xyy*eps ) + 
					  e.yyx*map( pos + e.yyx*eps ) + 
					  e.yxy*map( pos + e.yxy*eps ) + 
					  e.xxx*map( pos + e.xxx*eps ) );
}

int AA = 3;

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    // camera movement	
	//iTime
	vec3 ro = vec3( 0/*cos(iTime)*/, 0, 2.0 /** sin(iTime)*/), ta = vec3(0), tot = vec3(0);
       
    for( int m=0; m<AA; m++ )
        for( int n=0; n<AA; n++ ) {
            // pixel coordinates
            vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5,
                p = (2.0*(fragCoord+o)-iResolution.xy)/iResolution.y;

            // create view ray
            vec3 rd = normalize(vec3(p.x, sin(iTime)*p.y, -0.5));

            // raymarch
            float tmax = 3.0, t = 0.0;
            
            for( int i=0; i<256; i++ ) {
                vec3 pos = ro + t*rd;
                float h = map(pos);
                
                if( h<1e-4 || t>tmax ) i=256;
                
                t += h;
            }

            // shading/lighting	
            vec3 col = vec3(0.0);
            
            if( t<tmax ) {
                vec3 pos = ro + t*rd, nor = calcNormal(pos);
                
                float dif = clamp( dot(nor,vec3(0.58)), 0.0, 1.0 ),
                    amb = (1.0 + dot(nor,vec3(0,1,0)))/2.0;
                
                col = vec3(0.2,0.3,0.4)*amb + vec3(0.8,0.7,0.5)*dif;
            }

            // gamma        
            tot += sqrt(col);    
        }
        
	fragColor = vec4( tot/float(AA*AA), 1.0 );
}
