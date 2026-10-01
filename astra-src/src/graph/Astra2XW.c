/*The file includes the C-functions:
     createGC         Change_Color     Open_Screen  close_Screen
     Open_Window      MoveArrow        PutColorName initDefaultColors
     PSASetForeground PSADrawRectangle PSADrawLine  PSADrawLString
and interfaces for FORTRAN calls
    initvm    drawvm  rectvm  erasrw  setlin  redraw
    psopen    psclos  pcurso  rcurso  textvm  textbf  pscom   setcolor
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h> /* Used for sleep(unsigned int) */
#include <stdlib.h> /* Used for system(const char*) */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <X11/cursorfont.h>
#include <sys/times.h>
#include <stdint.h>

typedef int32_t INT_;

void initvm_(INT_*, INT_*, INT_*, INT_*, INT_*, char*, INT_*);
void redraw_();
void erasrw_();
void resizeWindow(unsigned int, unsigned int);
void textvm_(INT_*, INT_*, char*, INT_*);
void textbf_(INT_*, INT_*, char*, INT_*);
void textnb_(INT_*, INT_*, char*, INT_*);
void putString(GC, INT_*, INT_*, char*, INT_*);
void createpixmap_(INT_*);
void changeGCcolor(GC, INT_*);
void setcolor_(INT_*);
void pscom_(char*, INT_*);
void psopen_(char*, INT_*, INT_*);
void psclose_();
void close_screen_();
void drawvm_(INT_*, INT_*, INT_*, INT_*, INT_*); 
void rectvm_(INT_*, INT_*, INT_*, INT_*, INT_*);
void cleare(INT_*, INT_*, INT_*, INT_*, INT_*); 
void cleare_(INT_*, INT_*, INT_*, INT_*, INT_*);
void GetRWgeometry(int*, int*);
void PSADrawRectangle(double, double, double, double);
void PSADrawLine(double, double, double, double);
void PSASetForeground(int);
void PSADrawLString(double, double, char*, int, double);

#define BORDER_WIDTH 2
#define NORMAL_WINDOW 0
#define POP_UP_WINDOW 1
#define DEFAULT_GEOMETRY NULL
#define DEFAULT_FONT1 "variable"
#define STDfont "8x13"
#define HGHfont "8x13bold"
#define F1sh 13
#define maxPixels 66
#define maxPixmaps 6
/* Pixmaps' allocation:
   Pixmaps[0]=theRootWindow  Pixmaps[3] ECR footprint
   Pixmaps[1] wall, etc.  Pixmaps[4] free
   Pixmaps[2] NBI footprint  
   Pixmaps[5] used for intermediate storage (normally empty)
*/

Display *theDisplay;
Window theRootWindow;
GC theGCA, hghGC, hgh_menuGC, hintGC;
Colormap theColormap;
Drawable Pixmaps[maxPixmaps] = {0, 0, 0, 0, 0, 0};
Drawable theIconPixmap;
Cursor thePauseCursor, theMenuCursor;
int Xmode=0;  /* Default: NoX - Batch mode */
int XWX, XWY, XWW, XWH, XCorrection=0, YCorrection=0;
int theScreen, theDepth, theWidth, theHeight, iconState;
int ColorNum, theCurrentColorNo;
int AstraColorNum[64] = 
{
     0, 0, /* Color 0,   White (no change allowed) */
     1, 0, /* Color 1,   Black (no change allowed) */
    50, 0, /* Color 2,   Curve 1  (default Red) */
     3, 0, /* Color 3,   Curve 2  (default Blue) */
    62, 0, /* Color 4,   Curve 3  (default VioletRed) */
    37, 0, /* Color 5,   Curve 4  (default MediumSeaGreen) */
     5, 0, /* Color 6,   Curve 5  (default Brown) */
    15, 0, /* Color 7,   Curve 6  (default DarkTurquoise) */
    20, 0, /* Color 8,   Curve 7  (default Goldenrod) */
    29, 0, /* Color 9,   Curve 8  (default LimeGreen) */
    64, 0, /* Color 10,  Curve 9  (default Yellow) */
    30, 0, /* Color 11,  Curve 10 (default DarkGreen) */
     7, 0, /* Color 12,  Curve 11 (default Coral) */
    48, 0, /* Color 13,  Curve 12 (default Pink) */
    30, 0, /* Color 14,  WarningColor (default Magenta) */
     1, 0, /* Color 15, */
     1, 0, /* Color 16, */
     1, 0, /* Color 17, */
     1, 0, /* Color 18, */
     1, 0, /* Color 19, */
     1, 0, /* Color 20, */
     1, 0, /* Color 21, */
     1, 0, /* Color 22, */
     1, 0, /* Color 23, */
     1, 0, /* Color 24, */
     1, 0, /* Color 25, */
     1, 0, /* Color 26, */
     1, 0, /* Color 27, */
     1, 0, /* Color 28, */
     1, 0, /* Color 29, */
    50, 0, /* Color 30  Message (default Red) */
    26, 0  /* Color 31  History (default LightBlue) */
 };
unsigned long theBlackPixel, theWhitePixel, theCurrentColor;
unsigned long thePixels[maxPixels]; 
char *theColorNames[maxPixels] = {
    "White",  "Black", "Aquamarine",  /* 0-2 */
    "Blue",  "BlueViolet", "Brown",   /* 3-5 */
    "CadetBlue", "Coral", "CornflowerBlue", /* 6-8 */
    "Cyan",  "DarkGreen", "DarkOliveGreen", /* 9-11*/
    "DarkOrchid", "DarkSlateBlue", "DarkSlateGrey", /*12-14*/
    "DarkTurquoise", "DimGrey", "Firebrick",  /*15-17*/
    "ForestGreen", "Gold",  "Goldenrod",  /*18-20*/
    "Grey",  "Green", "GreenYellow",  /*21-23*/
    "IndianRed", "Khaki", "LightBlue",  /*24-26*/
    "LightGrey", "LightSteelBlue", "LimeGreen",  /*27-29*/
    "Magenta", "Maroon", "MediumAquamarine", /*30-32*/
    "MediumBlue", "MediumForestGreen", "MediumGoldenrod", /*33-35*/
    "MediumOrchid", "MediumSeaGreen", "MediumSlateBlue", /*36-38*/
    "MediumSpringGreen", "MediumTurquoise", "MediumVioletRed", /*39-41*/
    "MidnightBlue", "Navy",  "Orange",  /*42-44*/
    "OrangeRed", "Orchid", "PaleGreen",  /*45-47*/
    "Pink",  "Plum",  "Red",   /*48-50*/
    "Salmon",  "SeaGreen", "Sienna",  /*51-53*/
    "SkyBlue", "SlateBlue", "SpringGreen",  /*54-56*/
    "SteelBlue", "Tan",  "Thistle",  /*57-59*/
    "Turquoise", "Violet", "VioletRed",  /*60-62*/
    "Wheat",  "Yellow", "YellowGreen",  /*63-65*/
};

double PS_xA= 72., PS_yA= 720., PSsc;
FILE *PSAfile;
/* FlagPSA = 0 (no hard copy), -1 (white), 1 (other colors) */
int FlagPSA=0, FigAcount=1; 

/********************************************************************/
void resizewindow_(unsigned int *width, unsigned int *height){
    printf("width=%d, height=%d\n", *width, *height);
    int result = XResizeWindow(theDisplay, theRootWindow, *width, *height);
    if (result == BadValue ) printf("   bad value!!!\n");
    if (result == BadWindow) printf("   bad window!!!\n");
    XClearArea(theDisplay, theRootWindow, 0, 0, 0, 0, True);
    XFlush(theDisplay);
}
		   
/********************************************************************/
void createGC(Window theWindow, GC *theNewGC, char *fontName){
    XGCValues theGCValues;
    XFontStruct *fontStruct; 
    unsigned long theValueMask;
    theValueMask = 0L;
    *theNewGC = XCreateGC(theDisplay, theWindow, theValueMask, &theGCValues);
    if (*theNewGC != 0){
        fontStruct = XLoadQueryFont(theDisplay, fontName);
        if (fontStruct != 0) XSetFont(theDisplay, *theNewGC, fontStruct->fid);
        XSetForeground(theDisplay, *theNewGC, theBlackPixel);
        XSetBackground(theDisplay, *theNewGC, theWhitePixel);
    }
}

/********************************************************************/
void Change_Color(GC anyGC, int frgColor, int bkgColor){
    XSetForeground(theDisplay, anyGC, thePixels[frgColor]);
    theCurrentColor = thePixels[frgColor];
    XSetBackground(theDisplay, anyGC, thePixels[bkgColor]);
}

/********************************************************************/
void Open_Screen(){
    theDisplay = XOpenDisplay(NULL);
    if (theDisplay == NULL){
        printf(">>> ERROR: Cannot establish a connection to the X Server %s\n", XDisplayName(NULL));
        return;
    }
    theScreen = DefaultScreen(theDisplay);
    theDepth = DefaultDepth(theDisplay, theScreen);
    iconState = 0;
    theBlackPixel = BlackPixel(theDisplay, theScreen);
    theWhitePixel = WhitePixel(theDisplay, theScreen);
    theColormap = XDefaultColormap(theDisplay, theScreen);
    thePauseCursor = XCreateFontCursor(theDisplay, 130);/*97 or 34*/
/*The display information:*/
    theWidth = DisplayWidth(theDisplay, theScreen);
    theHeight = DisplayHeight(theDisplay, theScreen);
}

/********************************************************************/
void close_screen_(){
    int i;
    for (i = 1; i < maxPixmaps-1; i++){
        if (Pixmaps[i]) XFreePixmap(theDisplay, Pixmaps[i]);
    }
    XDestroySubwindows(theDisplay, theRootWindow);
    XDestroyWindow(theDisplay, theRootWindow);
    XFreeCursor(theDisplay, thePauseCursor);
    XCloseDisplay(theDisplay);
}

/********************************************************************/
void makeIcon(Window theNewWindow){
  /* Get pre-defined two-color Icon Pixmap */

#define sqF1_width 21
#define sqF1_height 20
    static unsigned char sqF1_bits[] = {
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x01, 0xe0, 0xff,
        0x00, 0xe0, 0xff, 0xfc, 0xff, 0xff, 0xfc, 0xff, 0x7c, 0x3e,
        0xfc, 0x7c, 0x3e, 0xfc, 0x7d, 0x3e, 0xfc, 0x39, 0x3f, 0xfc,
        0x39, 0x3f, 0xfc, 0x3b, 0x3f, 0xfc, 0x93, 0x3f, 0xfc, 0x93,
        0x3f, 0xfc, 0x97, 0x3f, 0xfc, 0xc7, 0xff, 0xff, 0xc7, 0xff,
        0xff, 0xcf, 0x3f, 0xfc, 0xef, 0x3f, 0xfc, 0xff, 0xff, 0xff};

    theIconPixmap = XCreatePixmapFromBitmapData(theDisplay, theNewWindow, 
        sqF1_bits, sqF1_width, sqF1_height, //Use sqF1
        theWhitePixel, thePixels[50], theDepth); // Red on white

    return;
}

/********************************************************************/
Window Open_Window(int x, int y, int width, int height, int flag, char* theTitle, 
    int iconicState, Window theParent, Cursor theCursor){
  
    XSetWindowAttributes theWindowAttributes;
    XSizeHints theSizeHints;
    unsigned long theWindowMask;
    Window theNewWindow;
    XWMHints theWMHints;
    XClassHint theClassHint;

    theWindowAttributes.border_pixel = theBlackPixel;
    theWindowAttributes.background_pixel = theWhitePixel;
    theWindowAttributes.cursor  = theCursor;
    if (flag == POP_UP_WINDOW){
        theWindowAttributes.override_redirect = True;
        theWindowAttributes.save_under = True;
        theWindowMask = CWBackPixel | CWBorderPixel | CWCursor | CWOverrideRedirect;
    }
    else{
        theWindowAttributes.override_redirect = False;
        theWindowMask = CWBackPixel | CWBorderPixel | CWCursor;
    }
    theNewWindow = XCreateWindow(theDisplay, theParent, 
        x, y, width, height, BORDER_WIDTH, theDepth, InputOutput, 
    CopyFromParent, theWindowMask, &theWindowAttributes);
    theWMHints.input = True;
    if (iconicState == 0){
        theWMHints.initial_state = NormalState;
    }
    else{
        theWMHints.initial_state= IconicState;
    }
    theWMHints.flags = InputHint | StateHint;
    makeIcon(theNewWindow);
    theWMHints.icon_pixmap = theIconPixmap;
    theWMHints.flags  = theWMHints.flags | IconPixmapHint;
    XSetWMHints(theDisplay, theNewWindow, &theWMHints);
    XStoreName (theDisplay, theNewWindow, theTitle);
    theSizeHints.flags = USPosition | PSize | PMinSize | PMaxSize;
    theSizeHints.x = x;
    theSizeHints.y = y;
    theSizeHints.width = width;
    theSizeHints.height = height;
    theSizeHints.min_width = width;  /* Resizing */
    theSizeHints.min_height = height; /* forbidden */
    theSizeHints.max_width = width;
    theSizeHints.max_height = height;
    XSetNormalHints(theDisplay, theNewWindow, &theSizeHints);
    XMapWindow(theDisplay, theNewWindow);
    XFlush(theDisplay);
    if (flag == NORMAL_WINDOW) sleep(1);
    return theNewWindow;
}

/********************************************************************/
void DrawArrow(Window wind, int Xx, int Xy, int color){
    Change_Color(theGCA, color, 0);
    XDrawLine(theDisplay, wind, theGCA, Xx, Xy, Xx, Xy+5);
    XDrawLine(theDisplay, wind, theGCA, Xx+1, Xy+1, Xx+1, Xy+5);
    XDrawLine(theDisplay, wind, theGCA, Xx-1, Xy+1, Xx-1, Xy+5);
    XDrawLine(theDisplay, wind, theGCA, Xx+2, Xy+3, Xx+2, Xy+5);
    XDrawLine(theDisplay, wind, theGCA, Xx-2, Xy+3, Xx-2, Xy+5);
    XDrawLine(theDisplay, wind, theGCA, Xx+3, Xy+5, Xx+3, Xy+5);
    XDrawLine(theDisplay, wind, theGCA, Xx-3, Xy+5, Xx-3, Xy+5);
}

/********************************************************************/
void MoveArrow(Window wind, int fromx, int fromy, int tox, int toy){
    if (fromx > 0){
        DrawArrow(wind, fromx, fromy, 0); /* white to erase */
    }

    if (tox > 0){
        DrawArrow(wind, tox, toy, 3); /* 3 - blue */
    }
    Change_Color(theGCA, 1, 0);
}

/********************************************************************/
/* (  x >= 0,      y >= 0)  - upper left corner  of the rectangle
 (x+w <= XWW, y+h <= XWH) - lower right corner of the rectangle
 the function puts a mark (circle) of the color 1 (red) 
   with a centre at the position (*x, *y), 
   preliminary it saves an image under the mark, 
   finally, it puts the image at the same place */
void PutColorName(Window wind, int ix, int iy, int iw, int num){
    int len, i;
    len = strlen(theColorNames[num]);
    i = (iw - 7*len)/2;
    XClearArea(theDisplay, wind, ix, iy, iw, 16, False);
    XDrawImageString(theDisplay, wind, hintGC, ix+i, iy+13, 
    theColorNames[num], len);
}

/********************************************************************/
void initDefaultColors(){
    XColor theRGBColor, theHardwareColor;
    int theStatus, i;
    if (theDepth > 1){
        ColorNum = maxPixels;
        for (i=0; i < maxPixels; i++){
            theStatus = XLookupColor(theDisplay, theColormap, theColorNames[i],
                &theRGBColor, &theHardwareColor);
            if (theStatus != 0){
                theStatus = XAllocColor(theDisplay, theColormap, &theHardwareColor);
                if (theStatus != 0){
                    thePixels[i] = theHardwareColor.pixel;
                }
                else{
                    thePixels[i] = theBlackPixel;
                }
            }
        }
    }
    else{  /* Monochrome system */
        for (i=0; i < maxPixels; i++){
            if (strcmp("White", theColorNames[i]) == 0){
                thePixels[i] = theWhitePixel;
            }
            else{
                thePixels[i] = theBlackPixel;
            }
        }
    }
}

/********************************************************************/
void initvm_(INT_ *x, INT_ *y, INT_ *wid, INT_ *hei, INT_ *LineWidth, char* Title, INT_ *titlen){
// Use the unix command "bitmap" to create the data

    Window openWindow();
    int i, ix, iy;

    Xmode = 1;

    if (sizeof(int) != 4) printf(" >>> Warning >>> Incompatibility in color table\n");
    Open_Screen();
    initDefaultColors();
    XWX = *x;
    XWY = *y;
    XWW = *wid;
    XWH = *hei;
 /* Ignore *x, *y and put the window in the top right corner */
    XWX = theWidth - XWW - 2*(BORDER_WIDTH + 5);
    XWY = 5; 
    theRootWindow = Open_Window(XWX, XWY, XWW, XWH, 0, Title, 
        iconState, RootWindow(theDisplay, theScreen), theMenuCursor);
    GetRWgeometry(&ix, &iy);
/* Each bloody XWindow Manager tries to do it in its own unique manner */
    XCorrection = ix - XWX;
    YCorrection = iy - XWY;
    createGC(theRootWindow, &hghGC, HGHfont);
    createGC(theRootWindow, &theGCA, STDfont);
    createGC(theRootWindow, &hgh_menuGC, HGHfont);
    createGC(theRootWindow, &hintGC, "variable");
    XSetLineAttributes(theDisplay, hgh_menuGC, 2L, LineSolid, CapRound, JoinRound);
    XSetLineAttributes(theDisplay, theGCA, *LineWidth, LineSolid, CapRound, JoinRound);
    XDrawRectangle(theDisplay, theRootWindow, theGCA, 0L, 0L, *wid-1L, *hei-1L);
    XWarpPointer(theDisplay, None, theRootWindow, 0, 0, 0, 0, *wid/2, *hei-120);

    Pixmaps[0] = theRootWindow;
}

/********************************************************************/
void createpixmap_(INT_ *id){
    INT_ i=0, j, k;
    if (Pixmaps[*id] == 0){
        Pixmaps[*id] = XCreatePixmap(theDisplay, theRootWindow, XWW, XWH, theDepth);
        j = XWW - 1;
        k = XWH - 1;
        cleare(id, &i, &i, &j, &k);
    }
}

/********************************************************************/
void redraw_(){
    XFlush(theDisplay);
}

/********************************************************************/
void erasrw_(){
// Clear the whole ASTRA GUI window
    XClearWindow(theDisplay, theRootWindow);
}

/********************************************************************/
void cleare(INT_ *id, INT_ *ix, INT_ *iy, INT_ *iw, INT_ *ih){
    cleare_(id, ix, iy, iw, ih);
}

void cleare_(INT_ *id, INT_ *ix, INT_ *iy, INT_ *iw, INT_ *ih){
   if (Pixmaps[*id] == 0) return;
   if (*id == 0){
       XSetForeground(theDisplay, theGCA, theWhitePixel);
   }
   else{
       XSetForeground(theDisplay, theGCA, theBlackPixel);
   }
   XFillRectangle(theDisplay, Pixmaps[*id], theGCA, *ix, *iy, *iw, *ih);
   XSetForeground(theDisplay, theGCA, theCurrentColor);
}

/********************************************************************/
void rcurso_(){
    XUndefineCursor(theDisplay, theRootWindow);
}

/********************************************************************/
void pcurso_(){
// Enter "xfd -center -fn cursor" to get all cursors available
    XDefineCursor(theDisplay, theRootWindow, thePauseCursor);
}

/********************************************************************/
void drawvm_(INT_ *id, INT_ *x1, INT_ *y1, INT_ *x2, INT_ *y2){
// Draw line from (x1, y1) to (x2, y2)
    int Xx1, Xy1, Xx2, Xy2;
    double dx1, dy1, dx2, dy2;
    Xx1 = *x1 + 10;
    Xx2 = *x2 + 10;
    Xy1 = *y1 + 10;
    Xy2 = *y2 + 10;
    if (*id) XSetForeground(theDisplay, theGCA, ~theCurrentColor);
    XDrawLine(theDisplay, Pixmaps[*id], theGCA, Xx1, Xy1, Xx2, Xy2);
    if (*id) XSetForeground(theDisplay, theGCA, theCurrentColor);
    if (FlagPSA > 0){
        dx1 = Xx1*PSsc;
        dx2 = Xx2*PSsc;
        dy1 = Xy1*PSsc;
        dy2 = Xy2*PSsc;
        PSADrawLine(dx1, dy1, dx2, dy2);
    }
}

/********************************************************************/
void drawcurve_(int *id, int *n, double *X, double *Y){
    int j, Xx1, Xy1, Xx2, Xy2;
    double dx1, dy1, dx2, dy2;

    if (*id) XSetForeground(theDisplay, theGCA, ~theCurrentColor);
    for (j=0; j < *n-1; j++){
        Xx1 = X[j]   + 10;
        Xy1 = Y[j]   + 10;
        Xx2 = X[j+1] + 10;
        Xy2 = Y[j+1] + 10;
        XDrawLine(theDisplay, Pixmaps[*id], theGCA, Xx1, Xy1, Xx2, Xy2);
    }
    if (*id) XSetForeground(theDisplay, theGCA, theCurrentColor);
    if (FlagPSA > 0){
        for (j=0; j < *n-1; j++){
            dx1 = (X[j]   + 10.)*PSsc;
            dy1 = (Y[j]   + 10.)*PSsc;
            dx2 = (X[j+1] + 10.)*PSsc;
            dy2 = (Y[j+1] + 10.)*PSsc;
            PSADrawLine(dx1, dy1, dx2, dy2);
        }
    }
}

/********************************************************************/
void rectvm_(INT_ *id, INT_ *x, INT_ *y, INT_ *width, INT_ *height){
    int Xx, Xy;
    unsigned int Xwidth, Xheight;
    double x1, y1, w1, h1;
    Xx = *x;
    Xy = *y;
    Xwidth  = *width;
    Xheight = *height;
    if (*id) XSetForeground(theDisplay, theGCA, ~theCurrentColor);
    XDrawRectangle(theDisplay, Pixmaps[*id], theGCA, Xx, Xy, Xwidth, Xheight);
    if (*id) XSetForeground(theDisplay, theGCA, theCurrentColor);
    if (FlagPSA > 0){
        x1 = Xx*PSsc;
        y1 = Xy*PSsc; 
        w1 = Xwidth*PSsc;
        h1 = (Xheight - 80)*PSsc;
        PSADrawRectangle(x1, y1, w1, h1);
    }
}

/********************************************************************/
void setcolor_(INT_ *clnumb){
    changeGCcolor(theGCA, clnumb);
}

/********************************************************************/
void changeGCcolor(GC gc_in, INT_ *clnumb){
    int inum;
    if (*clnumb < 32){
        inum = 2*(*clnumb);
    }
    else{
        inum = 62;
    }
    Change_Color(gc_in, AstraColorNum[inum], AstraColorNum[inum+1]);
    theCurrentColorNo = *clnumb;
    if (FlagPSA) PSASetForeground(inum);
}

/********************************************************************/
void textbf_(INT_ *x, INT_ *y, char *str, INT_ *str_len){
    putString(hghGC, x, y, str, str_len); 
}

/********************************************************************/
void textvm_(INT_ *x, INT_ *y, char *str, INT_ *str_len){
    putString(theGCA, x, y, str, str_len);
}

/********************************************************************/
void putString(GC gc_in, INT_ *x, INT_ *y, char *str, INT_ *str_len){
    int Xx, Xy, Xstlen;
    double x1, y1, psth;
    Xx = *x + 10;
    Xy = *y + 10;
    Xstlen = *str_len;
    XDrawImageString(theDisplay, theRootWindow, gc_in, Xx, Xy, str, Xstlen);
    if (FlagPSA > 0){
        x1 = Xx*PSsc;
        y1 = Xy*PSsc;
        psth = F1sh*PSsc;
        PSADrawLString(x1, y1, str, Xstlen, psth);
    }
}

/********************************************************************/
void textnb_(INT_ *x, INT_ *y, char *str, INT_ *str_len){
    int Xx, Xy, Xstlen;
    double x1, y1, psth;
    Xx = *x + 10;
    Xy = *y + 10;
    Xstlen = *str_len;
    XSetForeground(theDisplay, theGCA, ~theCurrentColor);
    XDrawString(theDisplay, Pixmaps[2], theGCA, Xx, Xy, str, Xstlen);
    XSetForeground(theDisplay, theGCA, theCurrentColor);
    if (FlagPSA > 0){
        x1 = Xx*PSsc;
        y1 = Xy*PSsc;
        psth = F1sh*PSsc;
        PSADrawLString(x1, y1, str, Xstlen, psth);
    }
}

/********************************************************************/
void pscom_(char *str, INT_ *str_len){
    int j;
    if (FlagPSA > 0){
        fprintf(PSAfile, "%% ");
        for (j=0; j < *str_len; j++) fprintf(PSAfile, "%c", str[j]);
        fprintf(PSAfile, "\n");
    }
    return;
}

static int icount = 0;

/********************************************************************/
void psopen_(char *PSname, INT_ *sty, INT_ *iret){
  
    if (FlagPSA){
        *iret = -1;
        return;
    }

    PSAfile = fopen(PSname, "w");

    if (PSAfile != NULL){
        FlagPSA = 1;
        *iret = 0;
        fprintf(PSAfile, "%%!PS-Adobe-2.0\n");
        fprintf(PSAfile, "/Lshow {exch dup scale show dup scale}def\n");
        fprintf(PSAfile, "/Rshow {dup stringwidth pop 3 -1 roll\n");
        fprintf(PSAfile, "dup 3 -1 roll mul neg 0 rmoveto\n");
        fprintf(PSAfile, "dup scale show dup scale}def\n");
        fprintf(PSAfile, "/Mshow {dup stringwidth pop 3 -1 roll\n");
        fprintf(PSAfile, "dup 3 -1 roll -0.5 mul mul 0 rmoveto\n");
        fprintf(PSAfile, "dup scale show dup scale}def\n");
        fprintf(PSAfile, "%%Page: %d %d\n", icount, icount);

        switch(*sty){
        case 0:
            PSsc = 0.7;
            fprintf(PSAfile, "%%Figure %2d\n %g setlinewidth\n", FigAcount, 0.4);
            fprintf(PSAfile, "/Courier-Bold findfont 1.02 scalefont setfont\n");
            break;
        case 1:
            PSsc = 0.98;
            fprintf(PSAfile, "%%Figure %2d\n %g setlinewidth\n", FigAcount, 0.5);
            fprintf(PSAfile, "756 0 translate\n");
            fprintf(PSAfile, "90 rotate\n");
            fprintf(PSAfile, "/Courier-Bold findfont 1.02 scalefont setfont\n");
            break;
        default:
            printf(" >>> PSOPEN >>> Unknown style option: %d\n", *sty);
        }
        PSASetForeground(0);
        FigAcount++;
    }
    else{
        *iret = 1;
    }
    return;
}

/********************************************************************/
void psclose_(){
    if (FlagPSA){
        fprintf(PSAfile, "stroke\nshowpage\n\n");
        FlagPSA = 0;
	fclose(PSAfile);
    }
    return;
}

/********************************************************************/
void PSASetForeground(int i){
    XColor theRGBColor;
    float r, g, b;
    fprintf(PSAfile, "stroke\n");
    theRGBColor.pixel = thePixels[AstraColorNum[i]];
    XQueryColor(theDisplay, theColormap, &theRGBColor);
    if ((theRGBColor.red & theRGBColor.green & theRGBColor.blue) != 0xffff){
        FlagPSA = 1;
        r = theRGBColor.red;
        g = theRGBColor.green;
        b = theRGBColor.blue; 
        r = r/0xffff;
        g = g/0xffff;
        b = b/0xffff;
        fprintf(PSAfile, "%5.3f %5.3f %5.3f setrgbcolor\n", r, g, b);
    }
    else{
        FlagPSA = -1;
    }
}

/********************************************************************/
void PSADrawRectangle(double x, double y, double w, double h){
    fprintf(PSAfile, "newpath %%Rectangle\n");
    fprintf(PSAfile, "%10.3e %10.3e  moveto\n", x+PS_xA, PS_yA-y);
    fprintf(PSAfile, "%10.3e %10.3e %10.3e %10.3e rlineto rlineto\n",  w, 0., 0., -h);
    fprintf(PSAfile, "%10.3e %10.3e %10.3e %10.3e rlineto rlineto\n", -w, 0., 0.,  h);
    fprintf(PSAfile, "stroke\n");
}

/********************************************************************/
void PSADrawLine(double x1, double y1, double x2, double y2){
    fprintf(PSAfile, "%11.4e %11.4e %11.4e %11.4e moveto lineto\n",
        x2+PS_xA, PS_yA-y2, x1+PS_xA, PS_yA-y1);
}

/********************************************************************/
void PSADrawLString(double x, double y, char *ln, int n, double psfsc){
    char tlin[256];
    double psfsc1;
    int i;
    psfsc1 = 1./psfsc;
    for (i=0; i <= n-1 ; ++i){
        *(tlin+i) = *(ln+i);
        *(tlin+n) = '\0';
    }
    fprintf(PSAfile, "%10.3e %10.3e moveto\n", x+PS_xA, PS_yA-y);
    fprintf(PSAfile, "%15.8e %15.8e (%s) Lshow\n", psfsc1, psfsc, tlin);
}
