#define iabs(x) ((x)<0 ? -(x) : (x))
#define dabs(x) ((x)<0.0 ? -(x) : (x))

#define isgn(x) ((x)<0 ? -1 : 1)
#define sgn(x) ((x)<0.0 ? -1.0 : 1.0)
#define dsgn(x) ((x)<0.0 ? -1.0 : 1.0)

#define SIGN(a,b) ((b) >= 0.0 ? dabs(a) : -dabs(a))

#define step(x)   ((x)<0.0 ? 0.0 : (x))

#define min(x,y) ((x)<(y) ? (x) : (y))
#define max(x,y) ((x)>=(y) ? (x) : (y))

#define mod(i,n) ((i)>=0 ? (i)%(n) : ((n)+(i)%(n))%(n))
#define mod1(x)   ((x)<0 ? (x)+1.0 : (x)-((int) (x)))

#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

#define swap(a,b) tempr=(a); (a)=(b); (b)=tempr
