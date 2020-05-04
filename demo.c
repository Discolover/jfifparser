#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <GL/glut.h>

#define SOI  0xFFD8
#define APP0 0xFFE0
#define DHT  0xFFC4
#define SOF0 0xFFC0
#define SOS  0xFFDA
#define DQT  0xFFDB
#define EOI  0xFFD9

/*macros*/
#define byte() (LENGTH -= 1, *PTRTOBUFF++)
#define word() (LENGTH -= 2, (((*PTRTOBUFF++ << 8) | 0) | *PTRTOBUFF++))
#define length() (LENGTH = word(), LENGTH -= 2)
#define round(num) (((num) > 0) ? (int)((num) + 0.5) : (int)((num) - 0.5))
#define set(array, nbytes)                  \
    do {                                    \
        memcpy(array, PTRTOBUFF, nbytes);   \
        PTRTOBUFF += nbytes;                \
        LENGTH -= nbytes;                   \
    } while (0)
/* */

struct JfifData {
    char id[5];
    uint16_t version;
    uint8_t units;
    uint16_t xden, yden;
    uint8_t xthumb, ythumb;
    uint8_t *thumb;
} jfif;

struct FrameData {
    uint8_t precision;	
    uint16_t nlines;
    uint16_t nsamples;
    uint8_t ncomps;
    uint8_t compid; 
    uint8_t horzsampfact;
    uint8_t vertsampfact;
    uint8_t quanttblid;
} frame;

struct ScanData {
    uint8_t ncomps;	
    uint8_t compsel;
    uint8_t dc_ac;
    uint8_t s;
    uint8_t e;
    uint8_t h_l;
} scan;

/*prototypes*/
void buf_file(char *name);
int check_jfiffile();
void init_jfif();
void decode_image();
uint16_t get_marker();
void init_quanttbl();
void init_hufftbls();
void init_frame();
void init_scan();
void process_entropy();
void init_draw();
void init_du(int *du);
int get_magnitude(int nbits);
uint8_t get_symbol(int class);
uint8_t* try_symbol(int class, int len, uint16_t code);
int get_bit();
void unzz_du(int *du);
void deqnt_du(int *du);
void idct_du(int *du);
int idct_xy(int x, int y, int *du);
void shift_du(int *du);
void draw_du(int *du, int xpos, int ypos);
void die(const char *fmt, ...);
void draw_image();
void get_info();
/**/

uint8_t *BUFF, *PTRTOBUFF; 
uint16_t LENGTH;

uint8_t numbers[2][16];
uint16_t *codes[2][16];
uint8_t *symbols[2][16]; 
uint8_t qtbl[64];

/*zigzag order*/
int zz[64] = {
     0,  1,  5,  6, 14, 15, 27, 28, 
     2,  4,  7, 13, 16, 26, 29, 42, 
     3,  8, 12, 17, 25, 30, 41, 43,
     9, 11, 18, 24, 31, 40, 44, 53, 
    10, 19, 23, 32, 39, 45, 52, 54,
    20, 22, 33, 38, 46, 51, 55, 60,
    21, 34, 37, 47, 50, 56, 59, 61,
    35, 36, 48, 49, 57, 58, 62, 63
};

/*DEBUG FUNCTIONS*/
void
print_quanttbl()
{
	int i;	

	printf("Quantization table:\n");
	for (i = 0; i < 64; ++i) {
		printf("%3d ", qtbl[i]);
		if (!((i+1) % 8)) putchar('\n');
	}
}

void
print_hufftbl(int class)
{
	int i, j;
	uint8_t **syms, *nums;
	uint16_t **cds;

	nums = numbers[class];
	syms = symbols[class];
	cds = codes[class];

	printf("\nClass %s\n", (class) ? "AC" : "DC");
	printf("symbol\tcode\n");
	for (i = 0; i < 16; ++i, ++cds, ++syms) {
		if (nums[i]) printf("Length %d:\n", i + 1);
		for (j = 0; j < nums[i]; ++j) {
			printf("%X", (*syms)[j]);
			printf("\t%X\n", (*cds)[j]);
		}
	}
}

void
print_du(int *du)
{
	int i;

	for (i = 0; i < 64; ++i) {
		printf("%2d ", du[i]); 
		if (!((i+1) % 8)) putchar('\n');
	}
	putchar('\n');
}
/*END DEBUG FUNCTIONS*/

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);	
    glutInitDisplayMode(GLUT_RGBA|GLUT_SINGLE);
	buf_file(argv[1]);	
	PTRTOBUFF = BUFF;
	if (!check_jfiffile())
	    die("is not jfif compatible file");
	init_jfif();
	decode_image();
	printf("--- end of image ---\n");
	get_info();
	draw_image();

	//int *p;
	//int sample[64] = {
	//	-26, -3, -6, 2, 2, -1, 0, 0,
	//	0, -2, -4, 1, 1, 0, 0, 0,
	//	-3, 1, 5, -1, -1, 0, 0, 0,
	//	-3, 1, 2, -1, 0, 0, 0, 0,
	//	1, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0
	//};

	//int testqtbl[64] = {
	//	16, 11, 10, 16, 24, 40, 51, 61,
	//	12, 12, 14, 19, 26, 58, 60, 55,
	//	14, 13, 16, 24, 40, 57, 69, 56,
	//	14, 17, 22, 29, 51, 87, 80, 62,
	//	18, 22, 37, 56, 68, 109, 103, 77,
	//	24, 35, 55, 64, 81, 104, 113, 92,
	//	49, 64, 78, 87, 103, 121, 120, 101,
	//	72, 92, 95, 98, 112, 100, 103, 99
	//};
	//for (p = testqtbl; p < testqtbl + 64; ++p)
	//	qtbl[p - testqtbl] = *p;
	//deqnt_du(sample);
	///*DEBUG*/print_du(sample);
	//idct_du(sample);
	///*DEBUG*/print_du(sample);
	//shift_du(sample);
	///*DEBUG*/print_du(sample);

	return 0;
}

void
buf_file(char *name)
{
    FILE *fp;			
    long flen;

    fp = fopen(name, "rb");
    fseek(fp, 0, SEEK_END);	
    flen = ftell(fp);
    rewind(fp);

    BUFF = malloc((flen + 1) * sizeof *BUFF);
    fread(BUFF, flen, 1, fp);
    fclose(fp);
}

int
check_jfiffile()
{
    return word() == SOI && word() == APP0;
}

void 
init_jfif()
{
   uint8_t n;

   length(); 
   set(jfif.id, sizeof jfif.id);
   jfif.version = word();
   jfif.units = byte();
   jfif.xden = word();
   jfif.yden = word();
   jfif.xthumb = byte();
   jfif.ythumb = byte();
   n = 3 * jfif.xthumb * jfif.ythumb;
   jfif.thumb = malloc(n * sizeof *jfif.thumb);
   set(jfif.thumb, n); 
}

void
decode_image()
{
    uint16_t m;

	printf("decode\n");

    while ((m = get_marker()) != EOI) {
        switch(m) {
        case DQT:
            init_quanttbl();
			print_quanttbl();
            break;
        case DHT:
            init_hufftbls();
            break;
        case SOF0:
            init_frame();
            break;
        case SOS:
            init_scan();
            break;
        default:
            die("Unknown marker %X", m);
        }
    }
}

uint16_t 
get_marker()
{
    uint8_t byte;	

    while ((byte = byte()) == 0xFF) 
        ;

	printf("> %X <\n", byte);

    return 0xFF00 | byte;
}

void
init_quanttbl()
{
    length();
    if (byte()) die("this program can have only 1 htbl OR this program support"
                        "only 8 bit precision");
    set(qtbl, sizeof qtbl);
}

void
init_hufftbls()
{
    uint8_t c_i, *nums, **syms, *p;
    uint16_t code, **cds;
    int i;

    length();	
    while (LENGTH) {
        c_i = byte();

        cds = codes[c_i >> 4];
        syms = symbols[c_i >> 4];
        nums = numbers[c_i >> 4];

        set(nums, 16);

        for (p = nums, code = 0; p < nums + 16; ++p, ++cds, ++syms) {
            *cds = malloc(*p * sizeof **cds);
            *syms = malloc(*p * sizeof **syms); 
            for (i = 0; i < *p; ++i, ++code) {
                (*cds)[i] = code;
                (*syms)[i] = byte();
            }
            code <<= 1;
        }
    }
	print_hufftbl(c_i >> 4);
}

void 
init_frame()
{
	uint8_t b;

    length();
    frame.precision = byte();
    frame.nlines = word();
    frame.nsamples = word();
    frame.ncomps = byte();
    frame.compid = byte(); 
	b = byte();
    frame.horzsampfact = b >> 4;
    frame.vertsampfact = b & 0x0f;
    frame.quanttblid = byte();
}

void
init_scan()
{
    length();
    set(&scan, sizeof scan);

    process_entropy();
}

void render()
{
}

/*ENTROPY PROCESSING*/
void
process_entropy()
{
    int j;
	int xdu, ydu, xndu, yndu;
    int *du;

	xndu = frame.nsamples / 8;
	yndu = frame.nlines / 8;
	init_draw();
    du = malloc(64 * sizeof *du);
	for (ydu = 0; ydu < yndu; ++ydu) {
		for (xdu = 0; xdu < xndu; ++xdu) {
			memset(du, 0, 64 * sizeof(*du));
			init_du(du);
			/*DEBUG*/print_du(du);
			unzz_du(du);
			/*DEBUG*/print_du(du);
			deqnt_du(du);
			/*DEBUG*/print_du(du);
			idct_du(du);
			/*DEBUG*/print_du(du);
			shift_du(du);
			/*DEBUG*/print_du(du);
			draw_du(du, xdu * 8, ydu * 8);
		}
    }
}

void
init_draw()
{
	glutInitWindowSize(frame.nsamples, frame.nlines);
    glutInitWindowPosition(0,0);
    glutCreateWindow("jpgviewer");
    glClearColor(1,1,1,0);
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, frame.nsamples, frame.nlines, 0.0);
}

void
init_du(int *du)
{
    int *p;
    uint8_t sym;
    static int pred;

    p = du;

    *p = get_magnitude(get_symbol(0)) + pred;
    pred = *p;

    for (p += 1; p < du + 64; ++p) {
        sym = get_symbol(1);

        if (sym == 0x00) break;
        //if (sym == 0xF0) {
        //    p += 0xF;
        //    continue;
        //}

        p += sym >> 4;
        *p = get_magnitude(sym & 0x0f);
    }
}

int
get_magnitude(int nbits)
{
    int val, plus, mask, i;

    if (nbits < 1) return 0;

	for (i = 0, mask = 0; i < nbits; ++i) {
		mask <<= 1;
		mask |= 1;
	}

    val = plus = get_bit();
    nbits -= 1;

    while (nbits) {
        val <<= 1;
        val |= get_bit();
        nbits -= 1;
    }

	if (!plus) val = ~val;

	val &= mask;
    return (plus) ? val : -val;
}

uint8_t
get_symbol(int class)
{
    int len;
    uint16_t code;
    uint8_t *psym;

    code = 0;
    len = 0;
    do {
        code <<= 1;
        code |= get_bit();
    } while (!(psym = try_symbol(class, len++, code)));

	printf("len: %d\n", len + 1); 
	printf("code: %X\n", code); 
	printf("symbol: %X\n", *psym); 

    return *psym;
}

uint8_t *
try_symbol(int class, int len, uint16_t code)
{
    uint16_t *cds, *p;
    uint8_t num;


    cds = codes[class][len];
    num = numbers[class][len];

    for (p = cds; p != NULL && p < cds + num; ++p)
        if (*p == code)
            return symbols[class][len] + (p - cds);		

    return NULL;
}

int
get_bit()
{
    static int offset = -1;
    static uint8_t byte;

    if (offset < 0) {
        offset = 7;
        byte = byte();
        /*discard 0x00 after 0xFF*/
        if (byte == 0xFF && byte() != 0x00) {
            /*probably this is a RESTART MARKER*/
            die("%X marker in entropy data", *(--PTRTOBUFF));
        }
    }

    return (byte >> offset--) & 1;
}
    
void
unzz_du(int *du)
{
    int cp[64], i;

    memcpy(cp, du, 64 * sizeof(*cp));

    for (i = 0; i < 64; ++i) {
        du[i] = cp[zz[i]];
    }
}
    
void
deqnt_du(int *du)
{
    uint8_t *p;

    for (p = qtbl; p < qtbl + 64; ++p, ++du)
        *du *= *p;
}

void
idct_du(int *du)
{
    int x, y;
    int cp[64];

    memcpy(cp, du, sizeof cp);

    for (x = 0; x < 8; ++x) {
        for (y = 0; y < 8; ++y) {
           *du++ = idct_xy(x, y, cp);
        }
    }
}

int
idct_xy(int x, int y, int *du)
{
    int u, v;
    double val, sum1, sum2;

    for (u = 0, sum1 = 0; u < 8; ++u) {
        for (v = 0, sum2 = 0; v < 8; ++v) {
            val = (u) ? 1 : (1 / sqrt(2));	
            val *= (v) ? 1 : (1 / sqrt(2));	
            val *= du[u * 8 + v];
            val *= cos(((2 * x + 1) * u * M_PI) / 16);
            val *= cos(((2 * y + 1) * v * M_PI) / 16);
            sum2 += val;
        }
        sum1 += sum2;
    }

    return round(0.25 * sum1);
}
    
void
shift_du(int *du)
{
    int *p;

    for (p = du; p < du + 64; ++p) {
        *p += 128;
        if (*p > 255) *p = 255;
    }
}

void
draw_du(int *du, int xpos, int ypos)
{
	int x, y;
	double r, g, b;
	glBegin(GL_POINTS);
	for (y = ypos; y < ypos + 8; ++y) {
		for (x = xpos; x < xpos + 8; ++x) {
			r = (double)(*du)/255; 
			g = (double)(*du)/255; 
			b = (double)(*du)/255; 
			glColor3f(r, g, b);
			glVertex2i(x, y);
			++du;
		}
	}
	glEnd();
}
/*END ENTROPY PROCESSING*/

/*DRAW IMAGE*/
void
draw_image()
{
	glFlush();
	glutDisplayFunc(render);
    glutMainLoop();
}
/*END DRAW IMAGE*/

/*GET INFO*/
void get_info()
{
    printf("units: %d\n", jfif.units);
    printf("xden: %d\n", jfif.xden);
    printf("yden: %d\n", jfif.yden);
    printf("nlines: %d\n", frame.nlines);
    printf("nsamples: %d\n", frame.nsamples);
    printf("ncomps: %d\n", frame.ncomps);
    printf("compid: %d\n", frame.compid); 
    printf("horzsampfact: %d\n", frame.horzsampfact);
    printf("vertsampfact: %d\n", frame.vertsampfact);
    printf("quanttblid: %d\n", frame.quanttblid);
    printf("ncomps: %d\n", scan.ncomps);	
    printf("compsel: %d\n", scan.compsel);
    printf("dc_ac: %d\n", scan.dc_ac);
    printf("s: %d\n", scan.s);
    printf("e: %d\n", scan.e);
    printf("h_l: %d\n", scan.h_l);
}
/*END GET INFO*/

void 
die(const char *fmt, ...)
{
    va_list args;

    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    fputc('\n', stderr);
    exit(1);
}
