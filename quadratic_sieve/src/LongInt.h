#ifndef __LONGINT_MYLIB__H__
#define __LONGINT_MYLIB__H__

#include <cstdarg>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

template<int n> class LNLgsB8n;

typedef LNLgsB8n<120> LongInt;

// ------------------------------------------------------------------------
// Fast! Binary numbers of 8*n bits maximum (no check on numbers limitation)
// ------------------------------------------------------------------------
template<int n> class LNLgsB8n_structures
{ public:
    LNLgsB8n<n> *_div, *_mod;
	LNLgsB8n_structures()
	{ _div=new LNLgsB8n<n>;
	  _mod=new LNLgsB8n<n>;
	}
	~LNLgsB8n_structures()
	{ delete _div;
	  delete _mod;
	}
};

template<int n> class LNLgsB8n
{ protected:
    static int i, j, k, l;
  public: 
	static LNLgsB8n_structures<n> SD; // static variables to save memory ant time
	~LNLgsB8n<n>() {}

    long r[n]; // the number in CC-256
    int sz;    // size of the vector with valid digits. Note: All algorithms doesnot check
	           // the value sz!!! If the size of the number is >sz, this will cause memory 
	           // problems
    int sign;  // 0+ 1- sign

	// - Constructors -------------------------------------------
	LNLgsB8n(void) { for(sz=n; sz; r[--sz]=0); sign=0; }
	
    LNLgsB8n(long t) { for(sz=n; sz; r[--sz]=0); sign=0; *this=t; }

    // - Initializators -----------------------------------------
	// LongInt = long;
    inline LNLgsB8n<n> & operator = (long value)
    { memset(r,0,sizeof(long)*(++sz));
	  sign=value>>31;
      r[sz=0]=value&0x7FFFFFFFL;
      Normalize();
      return *this;
    }

	// LongInt = LongInt;
    template<int size> LNLgsB8n<n> & operator = (LNLgsB8n<size> & b)
    { for(i=b.sz+1; i<=sz; ++i) r[i]=0;
      for(i=0; i<=b.sz; ++i) r[i]=b.r[i];
      sign=b.sign;
      sz=b.sz;
      return *this;
    }

	// 10 - system number in digits  LongInt="859047049456459" - is possible now!
	void operator = (const char * digits)
	{ memset(this, 0 ,sizeof(LNLgsB8n<n>));
	  for(i=0; digits[i]; ++i)
	    if(digits[i]>'9' || digits[i]<'0') break;
	    else
		{ for(j=0; j<=sz; ++j) r[j]*=10;
	      r[0]+=int(digits[i]-'0');
		  Normalize();
	    }
	}

	// LongInt = LongInt<size>;
    template<int size> void Copy(LNLgsB8n<size> & b)
    { i=n<b.sz?n:b.sz;
      for(j=sz+1; j<=i; ++j) r[j]=0;
      for(j=0; j<=i; ++j) r[j]=b.r[j];
      sign=b.sign;
      sz=i;
    }

    // - Normalization Function ---------------------------------
    inline void Normalize(void)
    { int i;
      for(i=0; i<=sz || (r[i]>>8); ++i) 
        r[i+1]+=r[i]>>8, r[i]&=0xFF;
      while(i && !r[i]) --i;
      sz=i;
    }

    // - Simple operators increment/decrement--------------------
	// ++LongInt;
    inline void operator ++(void)
    { if(sign) 
        { --r[i=0];
          if(r[i]<0)
            if(sz) 
              { while(r[i]<0) r[i]=0xFF, --r[++i];
                while(sz && !r[sz]) --sz;                
              }
              else
              { r[i]=1;
                sign^=1;
              }
        }
        else
        { ++r[i=0];
          while(r[i]>>8) r[i]&=0xFF, ++r[++i];
          if(i>sz) sz=i;
        }
    }

	// --LongInt;
    inline void operator --(void)
    { sign^=1;
      ++(*this);
      sign^=1;
    }


    // ----------------------------------------------------------
    // Shift functions
	// LongInt*=2;
    inline LNLgsB8n & p2(void) // fast multiplication on 2
    { for(i=j=0; i<=sz; ++i) r[i]=(r[i]<<1)|j, j=r[i]>>8, r[i]&=0xFF;
      if(j) r[++sz]=j;
      return *this;
    }

	// LongInt<<=shift; // binary shift to the left
    inline void operator <<= (int shift)
    { k=shift&0x7;
      if(k)
      { for(i=j=0; i<=sz; ++i) r[i]=(r[i]<<k)|j, j=r[i]>>8, r[i]&=0xFF;
        if(j) r[++sz]=j;
      }
      k=shift>>3; // x/8 == x>>3 , since 8==2^3
      if(k&&(sz||r[0]))
      { for(i=sz, j=sz+k; i>=0; --i, --j) r[j]=r[i];
        while(j>=0) r[j--]=0;
        sz+=k;
      }
    }

	// LongInt>>=shift; // binary shift to the right
    inline void operator >>= (int shift)
    { k=shift>>3; // x/8 == x>>3 , since 8==2^3
      if(k)
      { for(i=0, j=sz-k; i<=j; r[i++]=r[k++]);
        for(j=i; j<=sz; r[j++]=0);
        sz=i;
        while(!r[sz]&&sz) --sz;
      }
      k=shift&0x7;
      if(k)
      { j=r[0]>>k;  
        k=8-k;
        for(i=1,l=0; i<=sz; ++i,++l) r[l]=(r[i]<<k)|j, j=r[l]>>8, r[l]&=0xFF;
        if(j) r[sz]=j;
        else { r[sz]=0; if(sz) --sz; }
      }
    }

    // ----------------------------------------------------------
    // Show() - function
    // 0 - 8n decimal values, separated by space
    // 1 - 8n hexadecimal values, separated by space
    // 2 - minimum hexadecimal values, separated by space
    // 3 - minimum hexadecimal values
    // 4 - minimum binary values
	// 5 - 10-CS decimal system number output
	// 6 - 10-CS number only, without any brackets!
    inline void Save(FILE * out=stdout, int reg=1, bool perenos=true)
    { if(reg!=6) fprintf(out,"[%c ", (sign?'-':'+'));
	  if(reg==6 && sign) fprintf(out,"-");
	  if(reg==5 || reg==6)
		{ LNLgsB8n<n> p(*this);
	      long dig[n*3], s=0;
		  while(p!=0)
		  { dig[s++]=p%10;
		    p/=10;
		  }
		  if(s)
		  for(i=s-1; i>=0; --i)
			  fprintf(out,"%ld",dig[i]);
		  else fprintf(out,"0");

            // old version
			/*const digsize=n*3;
			long dig[digsize], s=0, szd=0, msz=0;
			long mult[digsize];
			for(i=0; i<digsize; ++i) mult[i]=dig[i]=0;
			mult[0]=1;

			for(i=0; i<=sz; ++i)
			{   for(j=0; j<=msz; ++j) dig[j]+=mult[j]*r[i];
				for(j=0; j<=szd+5; ++j) dig[j+1]+=dig[j]/10, dig[j]%=10;
				for(szd+=5; !dig[szd] && szd; --szd);
				for(j=0; j<=msz; ++j) mult[j]<<=8;
				for(j=0; j<=msz+5; ++j) mult[j+1]+=mult[j]/10, mult[j]%=10;
				for(msz+=5; !mult[msz] && msz; --msz);
			}
			for(i=szd; i>=0; --i) fprintf(out,"%d",dig[i]);
			*/
	    }
	  else
      for(i=(!reg||reg==1)?n-1:sz; i>=0; --i)
        { if(reg==0) fprintf(out,"%d ",r[i]);
          if(reg==1) fprintf(out,"%X ",r[i]);
          if(reg==2) fprintf(out,"%02X ",r[i]);
          if(reg==3) fprintf(out,"%X" ,r[i]);
          if(reg==4) for(j=7; j>=0; --j) fprintf(out,"%d",(r[i]>>j)&1);
        }
      if(reg!=6) fprintf(out,"]");
      if(perenos) fprintf(out,"\n");
    }

	// Output the number in decimal representation ( by the default on the Screen)
	inline void DecOutput(char*comm=NULL, FILE * out=stdout)
	{ if(comm)
	   { printf("%s: ", comm);
	     Save(out, 5);
	   }
	  else
		Save(out, 6, false); 
	}

    // ----------------------------------------------------------
    // Degree function 2^ret + (delta)
    inline int degree(void)
    { for(j=0, i=r[sz]>>1; i; i>>=1, ++j);
      j+=sz<<3;
      return j;
    }

    // Precise degree 2^double. Returns the logarithm of the number in base of 2
    inline double ddegree(void)
    { if(!sz && !r[0]) return 1.;
      i=degree();
      double delta=i, x, y, z, eps=pow(2., -1000.);
      for(j=i-1; j>=0; --j)
        if((r[j>>3]>>(j&0x7))&1)
        { x=y=0.;
          for(k=j+1; k<=i; ++k)
            if((r[k>>3]>>(k&0x7))&1) x+=pow(2., k-j);
          x=1./x;
          for(k=1, z=x; k<64 && x>eps; ++k)
            y+=((k&1)?1.:-1.)*x/(double)k, x*=z;
          y/=log(2.);
          delta+=y;
        }
      return delta;
    }

    // ----------------------------------------------------------
    // Comparison functions
    inline bool Compare (const LNLgsB8n<n> & b, bool less, bool equal)
    { if(sign^b.sign)
        if(sign) return false^less;
        else return true^less;
      if(sz>b.sz) return true^sign^less;
      if(sz<b.sz) return false^sign^less;
      for(i=sz; i>=0; --i)
        { if(r[i]>b.r[i]) return true^sign^less;
          if(r[i]<b.r[i]) return false^sign^less;
        }
      return false^equal;
    }

    inline bool operator > (const LNLgsB8n<n> & b)
    { return Compare(b, false, false); }

    inline bool operator < (const LNLgsB8n<n> & b)
    { return Compare(b, true, false); }

    inline bool operator >= (const LNLgsB8n<n> & b)
    { return Compare(b, false, true); }

    inline bool operator <= (const LNLgsB8n<n> & b)
    { return Compare(b, true, true); }

    inline bool operator != (const LNLgsB8n<n> & b)
    { return !operator==(b); }

    inline bool operator == (const LNLgsB8n<n> & b)
    { if(sign!=b.sign) return false;
      if(sz!=b.sz) return false;
      for(i=sz; i>=0; --i)
        if(r[i]!=b.r[i]) return false;
      return true;
    }

    // - Plus (Fast) --------------------------------------------
    inline void operator +=(LNLgsB8n<n> & b)
    { if(sign^b.sign) // if (-a)+(+b) or (+a)+(-b)
        { b.sign=sign;
          *this-=b;
          b.sign^=1;
          return;
        }      
      for(i=0; i<=b.sz;)
        { j=(r[i]+=b.r[i])>>8;
          r[i]&=0xFF;
          r[++i]+=j;
        }
      if(r[i]) while((j=r[i]>>8)) r[i]&=0xFF, r[++i]+=j;
        else --i;
      if(sz<i) sz=i;
    }

    // - Minus (Fast) -------------------------------------------
    inline void operator -=(LNLgsB8n<n> & b)
    { if(sign^b.sign)
        { b.sign=sign;
          *this+=b;
          b.sign^=1;
          return;
        }
      for(i=0; i<=b.sz;)
        { j=((r[i]+=0x100-b.r[i])>>8)-1;
          r[i]&=0xFF;
          r[++i]+=j;
        }
      while(r[i]<0 && i<=sz) r[i]+=0x100, --r[++i];
      if(r[i]<0)
        { sz=i;
          sign^=1;
          while(i--) r[i]=0xFF-r[i];
          ++r[++i];
          while(r[i]&0x100) r[i]&=0xFF, ++r[++i];
          r[sz]=0;
        }
      while(sz && !r[sz]) --sz;
    }

    // - Multiplication (Fast) ----------------------------------
	// LongInt *=LongInt;
    inline void operator *=(LNLgsB8n<n> & b)		
    { for(i=sz+b.sz; i>sz; r[i--]=l)
        for(j=sz, k=i-j, l=0; !(j<0 || k>b.sz); --j, ++k) l+=r[j]*b.r[k];
      for(         ; i>=0; r[i--]=l)
        for(j=i , k=i-j, l=0; !(j<0 || k>b.sz); --j, ++k) l+=r[j]*b.r[k];
      for(i=0, sz+=b.sz+2; i<=sz; ++i)
        r[i+1]+=r[i]>>8, r[i]&=0xFF;
      while(!r[sz] && sz) --sz;
      sign^=b.sign;      
    }

    // - Multiplication (Fast) ----------------------------------
	// LongInt *=int;
    inline void operator *=(int b)
    { for(i=0; i<=sz; ++i) r[i]*=b;
      Normalize();
    }

    // - Operations with return ---------------------------------
/*  // multioperator
    inline LNLgsB8n<n> multioperator( void (LNLgsB8n<n>::*op) 
          (LNLgsB8n<n>&), LNLgsB8n<n>& arg )
    { LNLgsB8n<n> q(*this); ((&q)->*op)(arg); return q; }
*/

    inline LNLgsB8n<n> operator *(LNLgsB8n<n> & b)
    { LNLgsB8n<n> q(*this); q*=b; return q; }

    inline LNLgsB8n<n> operator +(LNLgsB8n<n> & b)
    { LNLgsB8n<n> q(*this); q+=b; return q; }

    inline LNLgsB8n<n> operator -(LNLgsB8n<n> & b)
    { LNLgsB8n<n> q(*this); q-=b; return q; }

	inline LNLgsB8n<n> operator /(LNLgsB8n<n> & b)
	{ LNLgsB8n<n> div, mod;
	  DivMod(b, div, mod);
	  return div;
	}

	inline LNLgsB8n<n> operator %(LNLgsB8n<n> & b)
	{ LNLgsB8n<n> div, mod;
	  DivMod(b, div, mod);
	  return mod;
	}

    inline LNLgsB8n<n> operator *(int b)
    { LNLgsB8n<n> Res;
      Res=*this;
      Res*=b;
      return Res;
    }

	// - Extended Euclidian algorithm (Not Optimized Yet) --------------
	// gcd(a,b) = d = x*a + y*b
    inline void Extended_Euclid(LNLgsB8n<n> &a, LNLgsB8n<n> &b, LNLgsB8n<n> &d, LNLgsB8n<n> &x, LNLgsB8n<n> &y)
	{ if(b==0)
		{ memcpy(&d, &a, sizeof(LNLgsB8n<n>));
	      memset(&x,  0, sizeof(LNLgsB8n<n>));
	      memset(&y,  0, sizeof(LNLgsB8n<n>));
		  x.r[0]=1;
		  return;
		}
	  LNLgsB8n<n> div, mod;
	  a.DivMod(b, div, mod);
      Extended_Euclid(b, mod, d, y, x);
	  div*=x;
	  y-=div;
	}

    // - GCD-Greater Common Divisor (Not Optimized Yet) -------------
	// return gcd(a, b)
	inline LNLgsB8n<n> gcd(LNLgsB8n<n> & a, LNLgsB8n<n> & b)
	{ if(!b.sz && !b.r[0]) return a;
	  else return gcd(b, a % b);
	}

	// - Division and Modulo (Fast) ---------------------------------
	// Fast recurse: GetFixRate -> max between [to..from] such that ra>=rb*GetFixRate
	inline long GetFixRate(long size, long *ra, long *rb, long from, long to)
	{ if(from==to) return from;
          register long i=0, s=0, t=(from+to+1)>>1;
	  for(; i<=size; ++i)
	    s-=ra[i]-t*rb[i], s=s>0?(1+((s-1)>>8)):0;
	  if(s) return GetFixRate(size, ra, rb, from, t-1);
	  return GetFixRate(size, ra, rb, t, to);
	}

	inline void DivMod(LNLgsB8n<n> & b, LNLgsB8n<n> & div, LNLgsB8n<n> & mod)
	{ memcpy(&mod, this, sizeof(LNLgsB8n<n>));
	  memset(&div, 0, sizeof(LNLgsB8n<n>));
	  for(register long dif=mod.sz - b.sz, t, i, s; dif>=0; --dif)
	  { if(mod.sz < b.sz + dif) continue;
	    /*if(mod.sz && b.sz)
		   t= ((mod.r[mod.sz]<<8)| mod.r[mod.sz-1]) / (((b.r[mod.sz-dif])<<8)|b.r[mod.sz-dif-1]);
		else t=mod.r[mod.sz] / b.r[mod.sz-dif];
		*/
		if((div.r[dif]=GetFixRate(mod.sz-dif, mod.r+dif, b.r, 0, 255 ))) //(t>2?t-2:0), (t<253?t+1:255)
		  if(!div.sz) div.sz=dif;
		for(i=0, s=0; i<=b.sz; ++i)
		{ mod.r[i+dif]-=div.r[dif]*b.r[i]+s;
		  if(mod.r[i+dif]<0) s=1+((-mod.r[i+dif]-1)>>8), mod.r[i+dif]+=s<<8;
		    else s=0;
		}
		mod.r[i+dif]-=s;
		while(!mod.r[mod.sz] && mod.sz) --mod.sz;
	  }
	}

	// LongInt /=LongInt; - division and truncation
	inline void operator/=(LNLgsB8n<n> & b)
	{ DivMod(b, *SD._div, *SD._mod);
	  memcpy(this, SD._div, sizeof(LNLgsB8n<n>));
	}

	// LongInt %=LongInt; - modulo
	inline void operator%=(LNLgsB8n<n> & b)
	{ DivMod(b, *SD._div, *SD._mod);
	  memcpy(this, SD._mod, sizeof(LNLgsB8n<n>));
	}

	// - Power ----------------------------------------------------
	// Power function, based on the Tailor expansion of EXP and LN
	inline LNLgsB8n<n> powfn(double pwr)
	{ double deg=pwr*log(2.)*ddegree();
	  double bits=log(deg)/log(2.);
	  double pb=pow(2., (double)(long)(bits+1));
	  LNLgsB8n<n+65> p, S(1),T(1);
	  register int i,k;
	  for(i=(int)bits+1; i>-256; --i, pb/=2., p<<=1)
		  if(deg>pb) deg-=pb, p.r[0]|=1;
      S<<=128; // Precision is 2^128
	  T=p; 
	  for(k=2; !(T==0); k++)
	    S+=T, T*=p, T>>=256, T/=k;
	  S>>=256;
	  LNLgsB8n<n> res;
	  res=S;
	  ++res;
	  return res;
	}

	// LongInt /=int; - division and truncation (Fast!)
	inline void operator /=(long t)
	{ register double q=0;
	  for(i=sz; i>=0; --i)
		{ q=q*256.+(double)r[i];
		  r[i]=(long)(q/(double)t);
		  q-=(double)r[i]*(double)t;
		}
	  while(!r[sz] && sz) --sz;
	}

	// LongInt %=int; - modulo (Fast!)
	inline long operator %(long t)
	{ register long q=0;
	  for(i=sz; i>=0; --i)
	  { q=(q*8)%t;
	    q=(q*8)%t;
	    q=(q*4+r[i])%t;
	  }
	  return q;
	}

    // - Division by t-long ------------------------------------------------------
	// (Fast!) Divide *this by t if remainder is 0. Valid only for t of size <=2^23
	// return true if division by t successed.
	inline bool Divide_If_Divisible(long t)
	{ memcpy(SD._mod, this, sizeof(LNLgsB8n<n>));
	  long pp;
	  for(i=sz, k=0; i>=0; --i)
	    { pp=SD._mod->r[i];
		  SD._mod->r[i] =(int)( ((double)k*256.+(double)pp) / (double)t );
		  k = (k*8)%t;
		  k = (k*8)%t;
		  k = (k*4+pp)%t;
	    }
	  if(!k)
	  { while(!SD._mod->r[SD._mod->sz] && SD._mod->sz) --SD._mod->sz;
	    memcpy(this, SD._mod, sizeof(LNLgsB8n<n>));
		return true;
	  }
	  return false;
	}

	// Generate a random number of size around 2^bits
	void rand(int bits)
	{ memset(this, 0, sizeof(LNLgsB8n<n>));
	  for(i=0; i<bits; ++i)
		  r[i>>3]|=(std::rand()%2)<<(i&0x7);
	  r[bits>>3]|=1<<(bits&7);
	  sz=1+(bits>>3);
	  while(!r[sz] && sz) --sz;
	}
};

template<int n> int LNLgsB8n<n>::i;
template<int n> int LNLgsB8n<n>::j;
template<int n> int LNLgsB8n<n>::k;
template<int n> int LNLgsB8n<n>::l;
template<int n> LNLgsB8n_structures<n> LNLgsB8n<n>::SD;

#endif
