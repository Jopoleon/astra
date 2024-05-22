#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <X11/Xutil.h>
#ifndef INT8
#define INT_ int
#else
#define INT_ long
#endif

void taskmenu_(INT_*);
void xaxis(INT_*);
void xaxis_(INT_*);
void mvcursor_(INT_*, INT_*, INT_*);
void ProcessRootWindowEvent(XEvent*);
INT_ pollevent_(INT_*);
INT_ waitevent_(INT_*, INT_*, INT_*);
int  Cursor_in_Box();
INT_ Root_window_event(INT_*, int);
void stcopy(char*, char*, int);
int num2str(double, char*, int);
void iroundA(char*, int*, int*, int*);

extern Display *theDisplay;
extern Window theRootWindow;
extern Cursor theMenuCursor;
extern GC theGCA, hghGC, hgh_menuGC, hintGC;
extern int theScreen, theDepth, iconState, ColorNum;
extern int XWX, XWY, XWW, XWH, XCorrection, YCorrection, theHeight;
int XRW, YRW;
extern int Xmode;
extern int AstraColorNum[];
static int iact=-1, ibcursor=-1;

#define POLL_EV_MASK ( ButtonPressMask | KeyPressMask | ExposureMask |\
        StructureNotifyMask | FocusChangeMask | EnterWindowMask | LeaveWindowMask )
#define WRITE_ XDrawImageString(theDisplay, theWindow,
#define MVPOINTER_ XWarpPointer(theDisplay, None, theWindow, 0, 0, 0, 0,
#define LINGCA_ XDrawLine(theDisplay, theWindow, theGCA,

#define BTyof 0     /*button text vertical offset */
#define BTxof 0     /*button text horizontal offset */
#define nbuttons_max 36

int n_buttons=28;  /*main Astra  menu table */
char *MAMT[nbuttons_max];
char MAMK[nbuttons_max];
double width_ratio, height_ratio;
#define Button struct BUTTON
Button {int xleft, yup, xright, ydown, length; char *label;};
Button MAMB[nbuttons_max];

struct {int xcur, ycur;} Kevent;

void Change_Color(GC, int, int);
void Menu_Column(int, int, int, char*[], Button[]);
void Put_button(Button, GC);
void MoveArrow(Window, int, int, int, int);
void changeGCcolor(GC, INT_*);
void PutColorName(Window, int, int, int, int);
void stcopy(char*, char*, int);
extern int isascii(int); // 0 if the character is not ASCII, nonzero if it is ASCII
extern int isprint(int); // check if a character passed as the argument is a printable character or not
extern int isalnum(int); // checks whether a character is alphabet or number
int nextevent(INT_*, INT_*, INT_*, Button[], char[]);
int menubox_(char[], INT_*, double*, char[], INT_*, INT_*, INT_*);
int nbibox_(char[], char[], char[], INT_*, INT_*, INT_*, INT_*);
int layoutbox_(char[], char[], char[], INT_*, INT_*, INT_*, INT_*);
int ufilebox_(INT_*, char[], char[], char[]);
int FindBoxNum(int, int, int, int, int, int);
int GetEsc(XKeyEvent);
int GetValue(XKeyEvent, char[], int, char[], int*);
int GetKey(XKeyEvent, char[], int*);
Window Open_Window(int, int, int, int, int, char[], int, Window, Cursor);

/**********************************************************************/
void stcopy(char *s1, char *s2, int n){
/* Copy n characters from string s2 to string s1 
   equivalent to { strncpy(s1, s2, n); s1[n] = '\0'; }  */
    int i;
    for (i=0; i<n; ++i){
        *(s1+i) = *(s2+i);
    }
    *(s1+i) = '\0';
}

/**********************************************************************/
int num2str(double x, char *str, int str_len){
/* 
Input:
   x       : number
   str_len : string length
Output:
   str : string of length str_len
*/

    const int max_len=12;
    int i, ne, is, ll, k, i5;
    char ch[30];

// Initialise output string
    for (i=0; i<str_len; i++){
        str[i] = ' ';
    }
    str[str_len] = '\0';
    if (x == 0.){
        str[str_len-1] = '0';
        return 0;
    }
    sprintf(ch, "%+.*e", max_len, x);
    i = sscanf((ch + max_len + 4), "%d", &ne);

    if (ch[0] == '-'){
        is = 1;
    }
    else{
        is = 0;
    }
    ll = str_len - is;
    if (ll < 1){
        str[str_len-1] = '*';
        return 1;
    }
    k = max_len + 2;

    ch[0] = ch[1];
    for (i=1; i<k; i++){
        ch[i] = ch[i+2];
	k--;
    }

    ch[k] = '\0';
    ne += (1 - k);
    i5 = 0;
    if (k > ll){
        ne = ne + k - ll;
        k = ll;
        iroundA(ch, &k, &ne, &i5);
    }

    while (k) {
        if (ne >= 0) {
            if (ne + k <= ll) {
                for (i=str_len-1; ne > 0; ne--, str[i--] = '0');
                for (k--; k >= 0; str[i--] = ch[k--]);
                if (is) str[i] = '-';
                return 0;
            }
            else{
                i = 2;
                if (ne >  9) i++;
                if (ne > 99) i++;
                if (k + i > ll){
                    k--;
                    if (k == 0){
                        str[str_len-1] = '*';
                        return 1;
                    }
                    ne++;
                    iroundA(ch, &k, &ne, &i5);
                }
                else{
                    i = str_len - i;
                    sprintf(str + i, "e%d", ne);
                    for (--k; k >= 0; str[--i]=ch[k--]);
                    if (is) str[--i] = '-';
                    return 0;
                }
            }
        } 
        else{
            if (ll > -ne){
                if (k - ll){
                    i = k + ne;
                    if (i >= 0){
                        i = str_len;
                        while (ne){
                            k--;
                            i--;
                            str[i] = ch[k];
                            ne++;
                        }
                        i--;
                        str[i] = '.';
                        while (k){
                            k--;
                            i--;
                            str[i] = ch[k];
                        }
                    }
                    else{
                        i = str_len;
                        while (k){
                            k--;
                            i--;
                            str[i] = ch[k];
                            ne++;
                        }
                        while (ne){
                            i--;
                            str[i] = '0';
                            ne++;
                        }
                        i--;
                        str[i] = '.';
                    }
                    if (is){
                        i--;
                        str[i] = '-';
                    }
                    return 0;
                }
                else{
                    k--;
                    ne++;
                    iroundA(ch, &k, &ne, &i5);
                }
            }
            else{
                i = 3;
                if (ne <  -9) i = 4;
                if (ne < -99) i = 5;
                if (ll >= k + i){
                    i = str_len - i;
                    sprintf(str + i, "e%d", ne);
                    while (k){
                        i--;
                        k--;
                        str[i] = ch[k];
                    }
                    if (is){
                        i--;
                        str[i] = '-';
                    }
                    return 0;
                }
                k--;
                if (k == 0){
                    str[str_len-1] = '*';
                    return 1;
                }
                ne++;
                iroundA(ch, &k, &ne, &i5);
            }
        }
    }
    return 0;
}

/**********************************************************************/
void iroundA(char *ch, int *k, int *ne, int *i5){
    int i;
    i = *k - 1;
    if (ch[*k] >= 53 + *i5){
        ch[i]++;
        *i5 = 0;
        if (ch[i] == 53) *i5 = 1;
    }
    while (ch[i] == 58){
        if (i > 0){
            ch[i] = 48;
            i--;
            ch[i]++;
            *i5 = 0;
            if (ch[i] == 53) *i5 = 1;
            (*k)--;
            (*ne)++;
        }
        else{
            ch[i] = 49;
            (*ne)++;
        }
    }
    while (ch[(*k) - 1] == 48){
        if ((*k) > 1){
            (*k)--;
            (*ne)++;
        }
    }
}

/**********************************************************************/
void xaxis(INT_ *modex){
    xaxis_(modex);
}

void xaxis_(INT_ *modex){ /* Change of X-coordinate for plots */
    switch(*modex){
    case 0:
    case 1:
        MAMT[0] = "16*f(a)";
        MAMT[1] = "8*f(a)";
        break;
    case 2:
        MAMT[0] = "16*f(rho)";
        MAMT[1] = "8*f(rho)";
        break;
    case 3:
        MAMT[0] = "16*f(psi)";
        MAMT[1] = "8*f(psi)";
        break;
    }
    return;
}

/**********************************************************************/
void taskmenu_(INT_ *modex){
/* Draw menu table in Astra interactive mode, at the bottom of the main graphic window */
    int Xx, Xy, i, ixx=8, iyy, ny, ratio;
    int Wx, Wy;
    unsigned int Ww, Wh, Wb, Wd;
    Window theRW;
    char *BUTEXT[nbuttons_max] = {
        "16*f(a)", "8*f(a)", "8*f(psi)", "2*f(a,t)",
        "2*f(R,t)", "8*f(t)", "Equil", "Layout",
        "Refresh", "Style", "Next", "Back",
        "Variables", "Constants", "Grids", "Get X-axis",
        "Save log", "U-files", "Land PS", "Port PS",
        "Write data", "Type model", "Type data", "Test",
        "Run", "Step", "Quit", "Help"};
    char BUTKEY[nbuttons_max] = {
        '1', '2', '3', '4',
        '5', '6', '8', 'M',
        'R', '.', 'N', 'B',
        'V', 'C', 'D', 'X',
        'I', 'U', 'Q', 'G',
        'F', 'L', 'T', 'S',
        '\015', '\040', '\057', 'H'};

    XGetGeometry(theDisplay,theRootWindow, &theRW, &Wx, &Wy, &Ww, &Wh, &Wb, &Wd);

    width_ratio  = (double)Ww/660.;
    height_ratio = (double)Wh/550.;
 
/* Draw separating lines between plots and menu */
    Change_Color(theGCA, AstraColorNum[2], AstraColorNum[3]);
    ny = XWH - 128.*height_ratio;
    XDrawLine(theDisplay, theRootWindow, theGCA, 0, ny, XWW-1, ny);
    ny = XWH - 112.*height_ratio;
    XDrawLine(theDisplay, theRootWindow, theGCA, 0, ny  , XWW-1, ny);
    XDrawLine(theDisplay, theRootWindow, theGCA, 0, ny+1, XWW-1, ny+1);
    for (i=0; i<n_buttons; i++){
        MAMT[i] = BUTEXT[i];
        MAMK[i] = BUTKEY[i];
    }
    xaxis(modex);
    iyy = XWH - 74*height_ratio;
    ratio = (double)(XWW*4)/(double)n_buttons;
    for (i=0; i<n_buttons/4; i++){
        Menu_Column(ixx + i*ratio, iyy, 4, MAMT+4*i, MAMB+4*i);
    }
    if (ibcursor >= 0) Put_button(MAMB[ibcursor], hgh_menuGC);
/* Menu titles */
    Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
    Xx = ixx + 5*width_ratio;
    Xy = iyy - 5*height_ratio;
    XDrawImageString(theDisplay, theRootWindow, hghGC, Xx, Xy, "Graphic mode", 12);
    Xx += (int)(188*width_ratio);
    XDrawImageString(theDisplay, theRootWindow, hghGC, Xx, Xy, "Select", 6);
    Xx += (int)( 94*width_ratio);
    XDrawImageString(theDisplay, theRootWindow, hghGC, Xx, Xy, "Control", 7);
    Xx += (int)( 94*width_ratio);
    XDrawImageString(theDisplay, theRootWindow, hghGC, Xx, Xy, "In/Out", 6);
    Xx += (int)(188*width_ratio);
    XDrawImageString(theDisplay, theRootWindow, hghGC, Xx, Xy, "Status", 6);
}

/**********************************************************************/
void Menu_Column(int xm, int ym, int nbutt, char *met[], Button butt[]){
    int i, lenm=10, font_sym_width, font_sym_height, xwid, vert_sep=4, hor_sep=1;
    double F1sw=7.7, F1sh=13.;      /*font 1 symbol width, height */

    font_sym_width  = width_ratio *F1sw + 1;
    font_sym_height = height_ratio*F1sh + 1;
    xwid = lenm*font_sym_width + hor_sep;
    for (i=0; i<nbutt; i++){
        butt[i].xleft  = xm;
        butt[i].xright = butt[i].xleft + xwid;
        butt[i].yup    = ym + i*(font_sym_height + vert_sep);
        butt[i].ydown  = butt[i].yup + font_sym_height;
        butt[i].length = strlen(met[i]);
        butt[i].label  = met[i];
        Put_button(butt[i], theGCA);
    }
}

/**********************************************************************/
void Put_button(Button but, GC aGC){
/* Button drawing */
    int xl, xr, yd, yu;
    xl = but.xleft;
    xr = but.xright;
    yd = but.ydown;
    yu = but.yup;
    Change_Color(aGC, 1, 0);

    XClearArea(theDisplay, theRootWindow, xl-2, yu-2, xr-xl+4, yd-yu+4, False);
    XDrawImageString(theDisplay, theRootWindow, aGC, xl+BTxof+2, yd-BTyof-3, but.label, but.length);

    XDrawLine(theDisplay, theRootWindow, aGC, xl+2, yd  , xr-2, yd  ); // horizontal down
    XDrawLine(theDisplay, theRootWindow, aGC, xl+2, yu  , xr-2, yu  );
    XDrawLine(theDisplay, theRootWindow, aGC, xl  , yd-2, xl  , yu+2);
    XDrawLine(theDisplay, theRootWindow, aGC, xr  , yd-2, xr  , yu+2);

    XDrawLine(theDisplay, theRootWindow, aGC, xl  , yd-2, xl+2, yd  ); // bottom left
    XDrawLine(theDisplay, theRootWindow, aGC, xl  , yu+2, xl+2, yu  ); // top left
    XDrawLine(theDisplay, theRootWindow, aGC, xr-2, yu  , xr  , yu+2); // top right
    XDrawLine(theDisplay, theRootWindow, aGC, xr  , yd-2, xr-2, yd  ); // bottom right
}

/**********************************************************************/
void AstraEvent(){
   INT_ key;
   if (!Xmode) return;
   Root_window_event(&key, 1);
   return;
}

/**********************************************************************/
INT_ pollevent_(INT_ *key){
    return Root_window_event(key, 0);
}

/**********************************************************************/
INT_ Root_window_event(INT_ *key, int ii){
/* Polling for events */
    XEvent theEvent;
    XKeyEvent theKeyEvent;
    XComposeStatus theComposeStatus;
    KeySym theKeySym;
    int i, count, mode=QueuedAfterReading;
    int theKeyBufferMaxLen=4;
    char ch, theKeyBuffer[5];

    *key = 0;
    XSelectInput( theDisplay, theRootWindow,
        (ButtonPressMask | KeyPressMask | EnterWindowMask | StructureNotifyMask |\
        FocusChangeMask | ExposureMask | LeaveWindowMask) );
    count = XEventsQueued(theDisplay, mode);
    if (ii && count != 0){   /* Call from other processes */
        XNextEvent(theDisplay, &theEvent);
        if (theEvent.xany.window == theRootWindow) ProcessRootWindowEvent(&theEvent);
        XFlush(theDisplay);
        return 0;
    }
    if (count == 0){
        XFlush(theDisplay);
        return 0;
    }
    XNextEvent(theDisplay, &theEvent);
/* See other values of case in file /usr/include/X11/X.h */
    switch(theEvent.type){
    case ButtonPress:
        Kevent.xcur = theEvent.xbutton.x;
        Kevent.ycur = theEvent.xbutton.y;
        i = Cursor_in_Box();
        if (i >= 0 && i <= n_buttons) *key = (int)MAMK[i];
        if (*key == '\033'){
            *key = 'c';
            return 1; /* <Ctrl>+C */
        }
        if (theDepth == 1 &&  *key == 'A') *key = 0;
        return 0;
    case KeyPress:
        Kevent.xcur = theEvent.xbutton.x;
        Kevent.ycur = theEvent.xbutton.y;
        theKeyEvent = theEvent.xkey;
        count = XLookupString(&theKeyEvent, theKeyBuffer,
            theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
        if (count > theKeyBufferMaxLen){
            printf("Keyboard error:  Bufferlength =%d\n",count);
            break;
        }
        if (theKeyEvent.state & Mod1Mask){
            if (theKeySym <= 127){
                ch = theKeySym;
                *key = ch;
                return 2;
            }
        }
        if (theKeyEvent.state & ControlMask){
            if (isascii((int)theKeySym)){
                ch = theKeySym;
                *key = ch;
                return 1;
            }
        }
        theKeyBuffer[count] = '\0';
        *key = theKeyBuffer[0];
        if (theDepth == 1 && (*key == 'A' || *key == 'a') ) *key = 0;
        return 0;
    case Expose:
        if (theEvent.xexpose.count == 0)
        *key = (int)'R';
        return 0;
    case UnmapNotify:
        return 65005;
    case FocusIn:
        return 65006;
    case EnterNotify:
        return 65006;
    default:
        break;
    }
    return 0;
}

/**********************************************************************/
void GetRWgeometry(int *XRW, int *YRW){
    Window theRW, theCW;
    int Wx, Wy, iX, iY;
    unsigned int Ww, Wh, Wb, Wd;

    XGetGeometry(theDisplay,theRootWindow, &theRW, &Wx, &Wy,  &Ww, &Wh, &Wb, &Wd);
    XTranslateCoordinates(theDisplay, theRootWindow, theRW, 0, 0, &iX, &iY, &theCW);
    *XRW = iX - Wx - XCorrection;
    *YRW = iY - Wy - YCorrection;

    return;
}

/**********************************************************************/
INT_ waitevent_(INT_ *theKey, INT_ *xCursor, INT_ *yCursor){
/* Waiting for events from Astra in WAIT mode & and from dialog windows */
    INT_ i;
    i = nextevent(theKey, xCursor, yCursor, MAMB, MAMK);
    while (i == 0 && *theKey == 0 ){
        i = nextevent(theKey, xCursor, yCursor, MAMB, MAMK);
    }
    return i;
}

/*********************** Waiting events *******************************/
int nextevent(INT_ *theKey, INT_ *xCursor, INT_ *yCursor, Button mbb[], char mbk[]){
/* ibcursor -     current box # or -1
   iact     - highlighted box # or -1 */
    XEvent theEvent;
    XKeyEvent theKeyEvent;
    XComposeStatus theComposeStatus;
    KeySym theKeySym;
    int length, jbox, longKey, theKeyBufferMaxLen=4;
    char theKeyBuffer[5];

    XSelectInput(theDisplay, theRootWindow,
        (POLL_EV_MASK | PointerMotionMask | ButtonMotionMask | ButtonReleaseMask) );
    XNextEvent(theDisplay, &theEvent);

    switch(theEvent.type){

    case MapNotify:
        *theKey = 0;
        return 65004;

    case UnmapNotify:
        *theKey = 0;
        return 65003;

    case KeyPress:
        Kevent.xcur = theEvent.xbutton.x;
        Kevent.ycur = theEvent.xbutton.y;
        *xCursor = Kevent.xcur;
        *yCursor = Kevent.ycur;
        theKeyEvent = theEvent.xkey;
        length = XLookupString (&theKeyEvent, theKeyBuffer,
                theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
        if (length > theKeyBufferMaxLen || length < 0){
            printf("String translation error:  Bufferlength =%d\n",length);
            break;
        }
        if (length == 0){
            if (theKeySym >= 65000){ /* A special key (as <Ctrl>, <Alt>, <F1>, Home, etc) is pressed */
                *theKey = theKeySym;
                switch(theKeySym){
                case XK_Left:  return 65361;
                case XK_Up:    return 65362;
                case XK_Right: return 65363;
                case XK_Down:  return 65364;
                default: return 0;
                }
            }
            else{
                *theKey = theKeySym;
                return 65000;
            }
        }
        if (theKeyEvent.state & Mod1Mask){
            if (theKeySym <= 127){
                *theKey = theKeySym;
                return 2;
            }
        }
        if (theKeyEvent.state & ControlMask){
            if (isascii((int)theKeySym)){
                *theKey =  theKeySym;
                return 1;
            }
        }
        longKey = theKeyBuffer[0];
        theKeyBuffer[1] = '\0';
        if (length > 0){
            if (theKeySym >= ' ' && theKeySym <= '~' ){ /* ASCII 32, 126  */
                *theKey = longKey;
                if (theDepth == 1 && (theKeySym == 'A' || theKeySym == 'a')) *theKey = 0;
            }
            else{
                *theKey = longKey;
                if (longKey == 8)  return 65361; /* BS equiv <- */
            }
            return 65000;
        }
        break;

    case ButtonPress:
        Kevent.xcur = theEvent.xbutton.x;
        Kevent.ycur = theEvent.xbutton.y;
        *xCursor = Kevent.xcur;
        *yCursor = Kevent.ycur;
        jbox = Cursor_in_Box();
        if (jbox >= 0) *theKey = (int)mbk[jbox];
        if (theDepth == 1 && *theKey == 'A') *theKey = 0;        /*No color table*/
        if (*theKey == 015) Put_button(mbb[jbox], theGCA);  /*Undo highlight*/
        if (*theKey != 0) jbox = -1;
        switch(jbox){
        case -1: return 65001;
        case 24: return 65363;  /*   >    */
        case 25: return 65502;  /*  End   */
        case 26: return 65361;  /*   <    */
        case 27: return 65496;  /*  Home  */
        case 28: return 65364;  /*   >>   */
        case 29: return 65504;  /*  PgDn  */
        case 30: return 65362;  /*   <<   */
        case 31: return 65498;  /*  PgUp  */
        default : return 65001;
        }

    case ButtonRelease:
        Kevent.xcur = theEvent.xbutton.x;
        Kevent.ycur = theEvent.xbutton.y;
        *xCursor = Kevent.xcur;
        *yCursor = Kevent.ycur;
        return 65002;

    case MotionNotify:
        Kevent.xcur = theEvent.xbutton.x;
        Kevent.ycur = theEvent.xbutton.y;
        *xCursor = Kevent.xcur;
        *yCursor = Kevent.ycur;
        ibcursor = Cursor_in_Box();
        if (iact == ibcursor) return 65000;   /* No changes */
        if (iact >= 0){                       /* Hgh_box -> std */
            Put_button(mbb[iact], theGCA);
            iact = -1;
        }
        if (ibcursor >= 0){                   /* Std_box -> hgh */
            Put_button(mbb[ibcursor], hgh_menuGC);
            iact = ibcursor;
        }
        *theKey = 0;
        return 65000;

    case Expose:
        if (theEvent.xexpose.count == 0){
            longKey = 'R';  /* ASCII 82 */
        }
        else{
            longKey = 0;
        }
        *theKey = longKey;
        return 65000;

    case ConfigureNotify:
        *theKey = 0;
        return 0;
    }
    return 0;
}

/**********************************************************************/
int Cursor_in_Box(){
/* Returns:  menu box number, -1 if no match;
   Makes use of global structures
   Button MAMB[NMRMB];    Event Kevent; */
    int i;
    for (i=0; i<n_buttons; i++){
        if (Kevent.xcur >= MAMB[i].xleft && Kevent.xcur <= MAMB[i].xright &&
            Kevent.ycur >= MAMB[i].yup   && Kevent.ycur <= MAMB[i].ydown){
            return i;
        }
    }
    return -1;
}

/**********************************************************************/
void mvcursor_(INT_ *key, INT_ *ix, INT_ *iy){
/* Move cursor by one pixel */
    switch(*key){
    case 361:
        *ix--; /* Left  */
        break;
    case 362:
        *iy--; /* Up    */
        break;
    case 363:
        *ix++; /* Right */
        break;
    case 364:
        *iy++; /* Down  */
        break;
    }
    XWarpPointer(theDisplay, None, theRootWindow, 0, 0, 0, 0, *ix, *iy);
}

/**********************************************************************/
int menubox_(char title[], INT_ *arr_size, double *array, char varNames[],
	     INT_* nameLength, INT_ *id, INT_ *editable){
    Window theWindow;
    XEvent theEvent;
    const int num_str_len=6, var_name_len=6, /* length of value and name */
        wsym=8, hsym=13,          /* symbol width and height */
        n_columns=4,              /* # of columns */
        xshif=5, yshif=3;         /* table corner */
    int UpLeftx, UpLefty,         /* window corner location */
        Width, Height,            /* window size */
        ix, iy, wbox, hbox,       /* parameter box corner and size*/
        nparam,                   /* # of parameters (array size) */
        spos=0, vpos,             /* edit & edit start positions */
        ixa0=0, iya0=0, ixa, iya, /* arrow (textcursor) position */
        name_len, jbox, iret, oldparam, ixold, iyold,
        lline, xButton, yButton, i, ii=-1, ind, ind1, ind2, esc_flag,
        ihelp, selalb=0, icol, dcol, irow, drow;
    double valn, param;
    char value[10], ovalue[10], stri[10], vsym='=',
        grep_str[128], var_name[10], legend[128];

    name_len = *nameLength;
    lline = var_name_len + num_str_len + 1;
    hbox = hsym + 3;
    wbox = wsym*(lline + 1);
    vpos = wsym*(var_name_len + 1);
    nparam = *arr_size;

    Width = 2*xshif + n_columns*wbox - wsym;
    Height= 2*yshif + ((nparam - 1)/n_columns + 2)*hbox;
    if (*editable == 0) jbox = 0;
    if (*id == 6){
        jbox = 0;
        selalb = -1;
    }
    else{
        jbox = 1;
        selalb = 1;
    }

    GetRWgeometry (&XRW, &YRW);
    UpLeftx = 2;
    i = XRW - Width - 8;
    UpLefty = YRW;
    if (i > UpLeftx) UpLeftx = i;
    i = YRW + Height + 30;
    if (i > theHeight) UpLefty = theHeight - Height - 30;
    theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0, title, 0,
        RootWindow(theDisplay, theScreen), theMenuCursor);
    XSelectInput (theDisplay, theWindow, POLL_EV_MASK);
    ihelp = 0;

// Create table
    Change_Color(theGCA, 1, 0);        /* white background, black foreground */
    for (i=0; i<nparam; i++){
        ix = wbox*(i%n_columns) + xshif;
        iy = hbox*(i/n_columns) + yshif;
        ind = i*name_len;
        WRITE_ theGCA, ix, hsym+iy, varNames+ind, name_len);
        WRITE_ theGCA, ix+vpos-wsym, hsym+iy, &vsym, 1);
        valn = *(array+i);
        num2str(valn, value, num_str_len);
        WRITE_ theGCA, ix+vpos, hsym+iy, value, num_str_len);
    }

    oldparam = nparam;
    ixold = ix;
    iyold = iy;
    for (i=0; i<num_str_len; i++) ovalue[i] = value[i];
    if (jbox == 1 && selalb != -1 && ii < 0) MVPOINTER_ xshif+40, yshif+13);
    Change_Color(theGCA, 50, 0);
    for (i=1; i<n_columns; i++){
	ind = xshif + i*wbox - wsym/2;
        LINGCA_ ind, 0, ind, Height-hbox);
    }
    i = Height - hbox;
    LINGCA_ 0, i, Width, i);
    i--;
    LINGCA_ 0, i, Width, i);
    WRITE_ hghGC, xshif, Height-3, "OK",2);
    if (jbox == 0 && *editable == 0) MVPOINTER_ xshif+5, Height-5);
    if (*id == 6) MVPOINTER_ xshif+Width-50, Height-5);
    i = n_columns*wbox/wsym - 5;
    ind = 50;
    if (i < 50) ind = i;
    strcpy(legend, "/ Esc - done;");
    if (*id <= 4){
        strcat(legend, " Button/Tab - select;");
        strcat(legend, " Return/Tab - enter;");
        strcat(legend, " ? - quick help"); }
    else{
        strcat(legend, "     Button/Tab - select;");
        strcat(legend, "     Return/Tab - enter;" );
    }
    Change_Color(hintGC, 42, 0);
    ind = strlen(legend);
    if (*editable == 1) WRITE_ hintGC, xshif+2*wsym, Height-3, legend, ind);
    Change_Color(hintGC, 1, 0);
    if (*editable == 0) WRITE_ theGCA, xshif+4*wsym, Height-3,
        "    Information table.   No changes permitted.     ", ind);
    if (*id == 6){
        if (selalb ==  0) WRITE_ hintGC, xshif+8+(i-6)*wsym, Height-3,
	    " (Select all)", 13);
        if (selalb == -1) WRITE_ hintGC, xshif+8+(i-6)*wsym, Height-3,
	    "(Unselect 0)", 12);
    }
    Change_Color(theGCA, 1, 0);

    ix = wbox*((jbox - 1)%n_columns) + xshif;
    iy = hbox*((jbox - 1)/n_columns) + yshif;
    ind1 = (oldparam - 1)*name_len;
    ind2 = (jbox - 1)*name_len;
    WRITE_ theGCA, ixold, hsym+iyold, varNames+ind1, name_len);
    Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
    WRITE_ hghGC, ix, hsym+iy, varNames+ind2, name_len);
    Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
    WRITE_ theGCA, ixold+vpos-wsym, hsym+iyold, &vsym, 1);
    Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
    WRITE_ hghGC, ix+vpos-wsym, hsym+iy, &vsym, 1);
    valn = *(array+jbox-1);
    num2str(valn, value, num_str_len);
    Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
    WRITE_ theGCA, ixold+vpos, hsym+iyold, ovalue, num_str_len);
    Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
    WRITE_ hghGC, ix+vpos, hsym+iy, value, num_str_len);
    Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
    spos = -1;
    for (i=0; i<num_str_len; i++) stri[i] = ' ';
    ixa = ix + vpos;
    iya = iy + hsym + 2;
    MoveArrow(theWindow, ixa0, iya0, ixa, iya);
    ixa0 = ixa;
    iya0 = iya;

    while(1){
        XFlush(theDisplay);
        XNextEvent(theDisplay, &theEvent);
        if (theEvent.xany.window == theRootWindow){
            ProcessRootWindowEvent(&theEvent);
            continue;
        }

	esc_flag = 0;
        if (theEvent.type == ButtonPress){ // mouse
            xButton = theEvent.xbutton.x;
            yButton = theEvent.xbutton.y;
            ind = FindBoxNum(xButton-xshif, yButton-Height+hbox+1, 4*wsym, hbox, 1, 1);
	    if (*editable == 1){
                if (jbox){
                    oldparam = jbox;
                    ixold = ix;
                    iyold = iy;
                    if (spos > 0){
                        sscanf(stri, "%6lf", &param);
			*(array+jbox-1) = param;
                        num2str(param, value, num_str_len);
                        spos = -1;
                    }
                    for (i=0; i<num_str_len; i++) ovalue[i] = value[i];
                }
                jbox = FindBoxNum(xButton-xshif+wsym/2, 
                    yButton-yshif, wbox, hbox, n_columns, nparam);
            }
            if (ind){
                iret = -1;
                jbox = oldparam;
	 	esc_flag = 1;
            }
            else{
                if (jbox != 0){
                    ix = wbox*((jbox - 1)%n_columns) + xshif;
                    iy = hbox*((jbox - 1)/n_columns) + yshif;
                    ind1 = (oldparam - 1)*name_len;
                    ind2 = (jbox - 1)*name_len;
                    WRITE_ theGCA, ixold, hsym+iyold, varNames+ind1, name_len);
                    Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
                    WRITE_ hghGC, ix, hsym+iy, varNames+ind2, name_len);
                    Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
                    WRITE_ theGCA, ixold+vpos-wsym, hsym+iyold, &vsym, 1);
                    Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
                    WRITE_ hghGC, ix+vpos-wsym, hsym+iy, &vsym, 1);
                    valn = *(array+jbox-1);
                    num2str(valn, value, num_str_len);
                    Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
                    WRITE_ theGCA, ixold+vpos, hsym+iyold, ovalue, num_str_len);
                    Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
                    WRITE_ hghGC, ix+vpos, hsym+iy, value, num_str_len);
                    Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
                    spos = -1;
                    for (i=0; i<num_str_len; i++) stri[i] = ' ';
                    ixa = ix + vpos;
                    iya = iy + hsym + 2;
                    MoveArrow(theWindow, ixa0, iya0, ixa, iya);
                    ixa0 = ixa;
                    iya0 = iya;
                }
                continue;
            }
        } // Mouse button

        if (esc_flag == 0){
            if (jbox == 0){
                if (GetEsc(theEvent.xkey)){
                    break;
                }
                else{
                    continue;
                }
            }
            iret = GetValue(theEvent.xkey, stri, num_str_len, "0123456789.-+eE", &spos);
        }

        if (spos >= 0){
            ixa = ix + vpos + spos*wsym;
            iya = iy + hsym + 2;
            MoveArrow(theWindow, ixa0, iya0, ixa, iya);
            ixa0 = ixa;
            iya0 = iya;
        }
        if (iret && spos > 0){
            sscanf(stri, "%8lf", &param);
            *(array+jbox-1) = param;
            num2str(param, value, num_str_len);
	    //		    printf("num2str 1: %s\n", value);
            spos = -1;
        }

        if (iret == -1) break;
        if (iret == -2){ // Question mark: document usage of a selected variable
            for (i=0; i<name_len; i++){
                ii = name_len-i;
                strncpy(var_name, varNames+(jbox-1)*name_len+i, ii);
                if (var_name[0] != ' ') break;
            }
            var_name[ii] = '\0';
            ii--;
            while (var_name[ii] == ' '){
                var_name[ii] = '\0';
                ii--;
            }
            ii++;

            switch(*id){
            case 1:        // Variables
                strcpy (grep_str, "grep -i \"");
                strncat(grep_str, var_name, ii);
                strcat (grep_str, " \" main/variables.txt");
                break;
            case 2:        // Constants, usage in equ file (missing)
                strcpy (grep_str,"grep -i -w ");
                strncat(grep_str, var_name, ii);
                strcat (grep_str, " tmp/model.tmp");
                break;
            case 3:        // Time, grid control
                strcpy (grep_str, "grep -i \"");
                strncat(grep_str, var_name, ii);
                strcat (grep_str," \" main/internal.txt");
                break;
            case 4:
                strcpy (grep_str, "grep -w ");
                strncat(grep_str, var_name, ii);
                strcat (grep_str," tmp/model.txt");
                break;
            case 5:
                printf("Curve presentation\n");
                break;
            case 6:
                printf("Mark time slices\n");
                break;
	    }
            if (*id < 4){
                ii = system(grep_str);
                if (ii == 128) printf("The constant %s is not assigned\n", var_name);
            }
	}

        Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
        if (spos < 0){
            spos = -1;
            for (i=0; i<num_str_len; i++) stri[i] = ' ';
            WRITE_ hghGC, ix+vpos, hsym+iy, value, num_str_len);
        }
        else{
            WRITE_ hghGC, ix+vpos, hsym+iy, stri, num_str_len);
        }
        Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);

        if (iret >= 2){
            oldparam = jbox;
            ixold = ix;
            iyold = iy;
            for (i=0; i<num_str_len; i++) ovalue[i] = value[i];
	    if (iret == 2){
                jbox++;
                if (jbox > nparam) jbox = 1;
	    }
	    else{
                dcol = iret%10 - 2;
                icol = (jbox - 1)%n_columns + dcol;
                if (icol < 0) icol = 0;
                if (icol >= n_columns) icol--;
                drow = iret/10 - 2;
                irow = (jbox - 1)/n_columns + drow;
                if (irow < 0) irow = 0;
                if (irow > (nparam - 1)/n_columns) irow--;
                jbox = irow*n_columns + icol + 1;
                if (jbox == oldparam && icol == n_columns-1 && drow != -1) jbox++;
                if (jbox == oldparam && icol == 0       && drow !=  1) jbox--;
                if (jbox < 1) jbox = 1;
                if (jbox > nparam) jbox = nparam;
            }
            ix = wbox*((jbox - 1)%n_columns) + xshif;
            iy = hbox*((jbox - 1)/n_columns) + yshif;
            ind1 = (oldparam - 1)*name_len;
            ind2 = (jbox - 1)*name_len;
            WRITE_ theGCA, ixold, hsym+iyold, varNames+ind1, name_len);
            Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
            WRITE_ hghGC, ix, hsym+iy, varNames+ind2, name_len);
            Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
            WRITE_ theGCA, ixold+vpos-wsym, hsym+iyold, &vsym, 1);
            Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
            WRITE_ hghGC, ix+vpos-wsym, hsym+iy, &vsym, 1);
            valn = *(array+jbox-1);
            num2str(valn, value, num_str_len);
            Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
            WRITE_ theGCA, ixold+vpos, hsym+iyold, ovalue, num_str_len);
            Change_Color(hghGC, AstraColorNum[8], AstraColorNum[9]);
            WRITE_ hghGC, ix+vpos, hsym+iy, value, num_str_len);
            Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
            spos = -1;
            for (i=0; i<num_str_len; i++) stri[i] = ' ';
            ixa = ix + vpos;
            iya = iy + hsym + 2;
            MoveArrow(theWindow, ixa0, iya0, ixa, iya);
            ixa0 = ixa;
            iya0 = iya;
        }
    }

    XWarpPointer(theDisplay, None, theRootWindow, 0, 0, 0, 0, Kevent.xcur, Kevent.ycur);
    XSetInputFocus(theDisplay, theRootWindow, RevertToPointerRoot, theEvent.xkey.time);
    XDestroyWindow(theDisplay, theWindow);
    ibcursor = -1;
    XFlush(theDisplay);
    return 0;
}


/**********************************************************************/
int layoutbox_(char title[], char template[], char array[], INT_ *len,
    INT_ *nrows, INT_ *ngroup, INT_ *add_sep_line)

/*
Called from ASXWIN and ASTWIN (file src/for/surv.f90) , invoked by key "M" from IFKEY

Input:
    title        - Title of the table
    template     - string defining a structure of the table and its 1st line
                   1st line does not appear if all non-'|' symbols are spaces
    array        - data (numbers or strings) for input and output
    len          - length of "array" element according to description
                   in the calling routine (80 in the example below)
    nrows        - number of rows in a table to be created
    ngroup       - if > 0 distance (in rows) between blue separating lines
    add_sep_line - if > 0 separates bottom of the table with a fat blue line

Returned value:
    -1  - Error (window was not created)
     0  - Normal exit (table updated)
    >0  - Temporary exit (the returned value goes to the calling routine)
Call from FORTRAN
    Example 1:
        integer	        nrow, len
        character*80    title, template, array(10)
        len = 80
        nrow = 5        ! must be <= 10
        title = "Example"//char(0)
        template = " field1 |field2|    | f4 |       "//char(0)
        array(1) = " CumBol &  10.1& 2  &    &-1.E+5 "
        array(2) = "                                 "	! Empty string
        array(3) = "                                 "	! Empty string
        array(4) = "                                 "	! Empty string
        array(5) = "                                 "	! Empty string
Symbols in the positions "|" will not appear in the table
4th parameter must be 80
5th parameter must be <= 10
    call layoutbox(title,template,array,80,3,0,0)

More examples in the file "src/for/surv.f90", subroutines ASXWIN ASTWIN
*/

{
    XEvent theEvent;
    Window theWindow;
    const int xshif=5, yshif=3,    /* table corner */
        wsym=8, hsym=13,           /* symbol width and height */
        ncol_max=20;	            /* max / actual # of columns */

    int ind, jbox_old, hbox,	    /* current/old box No. & height */
        irow0=0, wbox,              /* number of chars in box */
        Width, Height,              /* defaut window sizes */
        ixold, iyold, icold, ixa0, iya0, arrdim,
        j_col, nwid[20], nsta[20], ii=0,
        icol, isym, irow, ixa,      /* arrow (textcursor) position */
        spos=0,                     /* abs. and rel. position of symbol in table */
        lvalue,	                    /* 1 if box was changed, 0 otherwise */
        fill_flag;
    int UpLeftx, UpLefty,           /* window corner location */
        nend, jbox, iret, ix, ix1, iy, i, j;
    char stri[132];

    j = strlen(title);
    if (j > 132){
        printf("LAYOUTBOX >>> Title\n\"%s\"\n           is too long\n", title);
        return -1;
    }
    if (*len > 132){
        printf("LAYOUTBOX >>> Table \"%s\" . Input string is too long\n", title);
        return -1;
    }
    nend = strlen(template);
    if (nend > 132){
        printf("LAYOUTBOX >>> Table \"%s\" . Requested width is too large\n", title);
        return -1;
    }
    for (j=0; j<nend ; j++){
	if (template[j] != '|' && template[j] != ' ') irow0 = 1;
	break;
    }

/**** Split "template" in blocks ****/
    j_col = 0;
    nsta[0] = 0;
    for (j=0; j<nend; j++){
        if (template[j] == '|'){
            nwid[j_col] = j - nsta[j_col];
            nsta[j_col+1] = j + 1;
	    j_col++;
	}
        if (j_col > ncol_max){
            printf("%s %d\n", "LAYOUTBOX >>>  Too many input fields ", j_col);
            return -1;
        }
    }
    nwid[j_col] = nend - nsta[j_col];
    j_col++;
    arrdim = *len;
    ixa0 = jbox_old = 0;
    hbox = hsym + 3;
    Width  = 2*xshif + nend*wsym;
    Height = 2*yshif + (*nrows + 4 + irow0)*hbox;
    if (Kevent.xcur+Kevent.ycur < 10){
        Kevent.xcur = 330;
        Kevent.ycur = 430;
    }

    GetRWgeometry (&XRW,&YRW);
    UpLeftx = 2;
    j = XRW - Width - 8;
    UpLefty = YRW;
    if (j > UpLeftx) UpLeftx = j;
    j = YRW + Height + 30;
    if (j > theHeight) UpLefty = theHeight - Height - 30;
    theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0, title, 0,
	RootWindow(theDisplay,theScreen), theMenuCursor);
    MVPOINTER_ xshif+wsym*nwid[0]/2, yshif+irow0*hbox+hsym);
    XSelectInput(theDisplay, theWindow, POLL_EV_MASK);

/*** Create table ***/
    Change_Color(theGCA, 3, 0); /* white background, blue  foreground */
    if (irow0){
/*** Type template in blue ***/
        ix = xshif;
        for (j=0; j<j_col; j++){
            WRITE_ theGCA, ix, hsym, template+nsta[j], nwid[j]);
            if (j+1 == j_col) break;
            ix += wsym*(nwid[j]+1);
        }
/*** Draw blue horizontal lines ***/
        j = hbox;
        LINGCA_ 0, j+4, Width, j+4);
    }
    if (*ngroup > 0)
       for (j=(*ngroup+irow0)*hbox; j+4+yshif<Height-2*hbox; j += (*ngroup)*hbox) LINGCA_ 0, j+4, Width, j+4);
    if (*add_sep_line > 0){
        j = Height - (*add_sep_line + 2)*hbox - 3;
        LINGCA_ 0, j, Width, j);
        ++j;
        LINGCA_ 0, j, Width, j);
    }
/*** Draw blue vertical lines ***/
    ix = xshif - 0.5*wsym;
    for (j=0; j<j_col-1; j++){
        ix += wsym*(nwid[j] + 1);
        ix1 = ix;
        if (nwid[j]   == 0) ix1 -= 2;
        if (nwid[j+1] == 0) ix1 += 2;
        LINGCA_ ix1, 0, ix1, Height-2*hbox);
    }
/*** Draw bottom line & comments ***/
    WRITE_ theGCA, xshif, Height-3-hbox, "Button - select,  <Tab>, <Ret>, Arrows - move", 45);
    WRITE_ theGCA, xshif+4*wsym, Height-4, "/<ESC> - done", 13);
    Change_Color(theGCA, 50, 0);  /* white background, red foreground */
    j = Height - 2*hbox;
    LINGCA_ 0, j, Width, j);
    j++;
    LINGCA_ 0, j, Width, j);
    Change_Color(hghGC, 50, 0);
    WRITE_ hghGC,  xshif, Height-4, " OK ", 4);
    XDrawRectangle(theDisplay, theWindow, hghGC, 5L, Height-17L, 30L, 16L);
    Change_Color(theGCA, 1, 0);
    Change_Color(hghGC , 1, 0);      /* white background, black foreground */
/*** Draw input data ***/
    for (i=0; i<*nrows; i++){
        iy = yshif + hsym + hbox*(i + irow0);
        for (j=isym=0, ix=xshif; j<j_col; j++){
	    if (j > 0){
                isym = nsta[j];
                ix += wsym*(nwid[j-1]+1);
            }
            WRITE_ theGCA, ix, iy, array+i*arrdim+isym, nwid[j]);
	}
    }
    if (jbox_old > 0){  /* after Expose event */
        isym = 0;
        if (icol > 1)  isym = nsta[icol-1];
	ix = xshif + wsym*isym;
	iy = yshif + hsym + hbox*(irow + irow0 - 1);
        WRITE_ hghGC, ix, iy, stri, wbox);
    }
    if (ii == 0){    /* Select Upper Left box on entry */
        jbox = jbox_old = irow = icol = icold = ii = 1;
        ind = isym = lvalue = 0;
        wbox = nwid[0];
	ix = ixold = ixa = ixa0 = xshif + wsym*isym;
	iy = iyold = iya0 = yshif + hsym + hbox*irow0;
        WRITE_ hghGC, ix, iy, array, wbox);
        strncpy(stri, array, wbox);
    }
    XFlush(theDisplay);

/*** Table Control ***/
    while(1){
        XNextEvent(theDisplay, &theEvent);
        if (theEvent.xany.window == theRootWindow){
            ProcessRootWindowEvent(&theEvent);
            continue;
        }
        if (theEvent.type == ButtonPress){
            ix = theEvent.xbutton.x;
            iy = theEvent.xbutton.y;
            irow = iy - yshif - hbox*irow0 - 2;
            if (irow >= 0){
                irow = irow/hbox + 1;
                i = (ix - xshif)/wsym;
                if (jbox_old > 0) strncpy(array+ind, stri, wbox);
                if (irow > *nrows) break; /* No selection or Exit */
                for (j=0; j<j_col; ){   /* Selection is made */
                    i -= nwid[j] + 1;
                    icol = ++j;
                    if (i < 0) break;
                }
                fill_flag = 1;
            }
        }
        else{
            fill_flag = 0;
	    iret = GetKey(theEvent.xkey, stri, &spos);
            switch(iret){
            case -1:       /* <Esc> */
                strncpy(array+ind, stri, wbox);
                break;
            case 0:
                if (spos == 0){           /* 1-st entry in the box */
                    for (j=1; j<132; j++) stri[j] = ' ';
                    stri[wbox] = '\0';
                }
                spos++;
                if (spos > wbox) spos = wbox;
                lvalue = 1;
                WRITE_ hghGC, ix, iy, stri, wbox);
                ixa = ix + spos*wsym;
                MoveArrow(theWindow, ixa0, iya0, ixa, iy);
                break;
            case 1: case 2:  /* <Tab> or <Ret> */
                icol++;
                if (icol > j_col) icol = 1;
                if (icol == 1){
                    irow++;
                    if (irow > *nrows) irow = 1;
                }
                if (spos != 0) strncpy(array+ind, stri, wbox);
                fill_flag = 1;
                break;
            case 3:             /* <Del> or <BackSpace> */
                spos--;
                if (spos >= 0){
                    strncpy(stri+spos, stri+spos+1, wbox-spos);
                    stri[wbox-1] = ' ';
                }
                else{
                    spos = 0;
                    strncpy(stri, array+ind, wbox);
                }
                WRITE_ hghGC, ix, iy, stri, wbox);
                ixa = ix + spos*wsym;
                MoveArrow(theWindow, ixa0, iya0, ixa, iy);
                break;
            case 12:              /* case XK_Up: */
                irow--;
                if (irow == 0) irow = *nrows;
                if (spos != 0 ) strncpy(array+ind, stri, wbox);
                fill_flag = 1;
                break;
            case 21:              /* case XK_Left: */
                spos--;
                if (spos < 0){
                    spos = 0;
                    icol--;
                    if (icol == 0){
                        icol = j_col;
                        irow--;
                    }
                    if (irow == 0) icol = irow = 1;
                    if (lvalue != 0) strncpy(array+ind, stri, wbox);
                    fill_flag = 1;
                }
                else{
                    ixa = ix + spos*wsym;
                    MoveArrow(theWindow, ixa0, iya0, ixa, iy);
                }
                break;
            case 23:              /* case XK_Right: */
                if (spos == 0 && lvalue == 0) strncpy(stri, array+ind, wbox);
                spos++;
                if (spos > wbox){ // Move to next box on the right
                    icol++;
                    if (icol > j_col) icol = 1;
                    if (icol == 1){
                        irow++;
                        if (irow > *nrows) irow = 1;
                    }
                    if (spos != 0) strncpy(array+ind, stri, wbox);
                    fill_flag = 1;
                }
                else{
                    ixa = ix + spos*wsym;
                    lvalue = 1;
                    MoveArrow(theWindow, ixa0, iya0, ixa, iy);
                }
                break;
            case 32:              /* case XK_Down: */
                irow++;
                if (irow > *nrows) irow = 1;
                if (spos != 0) strncpy(array+ind, stri, wbox);
                fill_flag = 1;
                break;
            default:
                break;
            }
            if (iret == -1) break; // Esc -> destroy child window
        }

// Fill box
        if (fill_flag == 1){
            isym = spos = 0;
            if (icol > 1) isym = nsta[icol-1];
            ix = ixa = xshif + wsym*isym;
            iy = yshif + hsym + hbox*(irow + irow0 - 1);
            ind = (irow-1)*arrdim + isym;
            wbox = nwid[icol-1];
            jbox = (irow - 1)*j_col + icol;
            if (jbox_old != 0) WRITE_ theGCA, ixold, iyold, stri, nwid[icold-1]);
            strncpy(stri, array+ind, wbox);
            WRITE_ hghGC, ix, iy, stri, wbox);
            MoveArrow(theWindow, ixa0, iya0, ixa, iy);
            icold = icol;
            ixold = ix;
            iyold = iy;
            jbox_old = jbox;
            lvalue = 0;
	    XFlush(theDisplay);
        }
        ixa0 = ixa;
        iya0 = iy;
    } /* End table control */

    XDestroyWindow(theDisplay, theWindow);
    XFlush(theDisplay);
    ibcursor = -1;
    return 0;
}

/**********************************************************************/
void Call_ifkey(int j){
    if (j == 47 ||                 /* "/" Exit -> */
        j == 13 ||                 /* <Enter> */
        j == 32 ||                 /* <Space> */
        j == 37 ||                 /* "%" */
        (j >= 48 && j <= 57  ) ||  /* "0" -> "9" */
        j == 66 ||                 /* "B" */
        j == 67 ||                 /* "C" */
        j == 68 ||                 /* "D" */
        j == 72 ||                 /* "H" */
        j == 78 ||                 /* "N" */
        j == 79 ||                 /* "O" */
        j == 82 ||                 /* "R" */
        j == 83 ||                 /* "S" */
        j == 87 ||                 /* "W" */
        j == 88                    /* "X" */
        ) ifkey_(&j);
    return;
}

/**********************************************************************/
void ProcessRootWindowEvent(XEvent *theEvent){
/* The function processes selected events in theRootWindow,
   e.g. redrawing in case of exposure. */
    int j, k;
    XKeyEvent theKeyEvent;
    XComposeStatus theComposeStatus;
    KeySym theKeySym;
    char theKeyBuffer[5];

    switch( theEvent->type ){
    case ButtonPress:
        Kevent.xcur = theEvent->xbutton.x;
        Kevent.ycur = theEvent->xbutton.y;
        k = Cursor_in_Box();         /* Returns box number */
        if (k >= 0 && k <= n_buttons) j = (int)MAMK[k];
        Call_ifkey(j);
        j = 0;
        break;
    case KeyPress:
        theKeyEvent = theEvent->xkey;
        j = XLookupString (&theKeyEvent, theKeyBuffer, 4, &theKeySym, &theComposeStatus);
        if (j == 0) return;          /* <Ret>, <Alt>, etc. is pressed */
        if ( j != 1 ){
            printf("String translation error:  Bufferlength =%d\n", j);
            return;
        }
        j = theKeySym;
        if (j > 96 && j < 123 ) j -= 32;
        Call_ifkey(j);
        j = 0;
        break;
    case Expose:
        j = 82;
        ifkey_(&j);
        j = 0;
        break;
    }
    return;
}

/**********************************************************************/
int nbibox_ (char title[], char template[], char array[], INT_ *len,
	     INT_ * nrows, INT_ *ngroup, INT_ *morow){
/* The same as "menubox", but the 1st column
   		is drawn in blue and closed for access
Input:	title	- Title of the table
	template  string defining a structure of the table and its 1st line
		  1st line does not appear if all non-'|' symbols are spaces
	array	- data (numbers or strings) for input and output
	len	- length of "array" element according to description
			in the calling routine (80 in the example below)
	nrows	- number of rows in a table to be created
	ngroup	- if > 0 distance (in rows) between blue separating lines
	morow	- if > 0 separates bottom of the table with a fat blue line

Called only from src/nbi/nbinj.f
Just a placeholder, if needed use the abstract "menubox"
*/
}

/**********************************************************************/
KeySym GetKeySym(XKeyEvent theKeyEvent){
    XComposeStatus theComposeStatus;
    KeySym theKeySym;
    int theKeyBufferMaxLen = 4;
    char theKeyBuffer[5];
    XLookupString(&theKeyEvent, theKeyBuffer, theKeyBufferMaxLen,
        &theKeySym, &theComposeStatus);
    return theKeySym;
}

/**********************************************************************/
int GetEsc(XKeyEvent theKeyEvent){
    KeySym theKeySym = GetKeySym(theKeyEvent);
    return (theKeySym == XK_Escape);
}

/**********************************************************************/
int GetKey(XKeyEvent theKeyEvent, char str[], int *pos){
    KeySym theKeySym = GetKeySym(theKeyEvent);

    switch(theKeySym){
    case XK_BackSpace:
    case XK_Delete:    return 3;

    case XK_question:  return -2;
    case XK_Escape:    return -1;
    case XK_Return:
    case XK_KP_Enter:  return 1;
    case XK_Tab:       return 2;
    case XK_R7:        return 11;
    case XK_Left:      return 21;
    case XK_Up:        return 12;
    case XK_Right:     return 23;
    case XK_Down:      return 32;
    case XK_R9:        return 13;
    case XK_R13:       return 31;
    case XK_R15:       return 33;

    default:
        if (!isprint(theKeySym)){
            return 99;
        }
        else{
            str[*pos] = theKeySym;
            break;
        }
    }
    return 0;
}

/**********************************************************************/
int GetName(XKeyEvent theKeyEvent, char str[], int lvalue, int *pos){
    int i;
    char lsym;
    KeySym theKeySym = GetKeySym(theKeyEvent);

    switch(theKeySym){
    case XK_BackSpace:
    case XK_Delete:
        *pos -= 1;
        if (*pos >= 0) str[*pos] = ' ';
        break;

    case XK_question: return -2;
    case XK_Escape:   return -1;
    case XK_Return:
    case XK_KP_Enter: return 1;
    case XK_Tab:      return 2;
    case XK_R7:       return 11;
    case XK_Left:     return 21;
    case XK_Up:       return 12;
    case XK_Right:    return 23;
    case XK_Down:     return 32;
    case XK_R9:       return 13;
    case XK_R13:      return 31;
    case XK_R15:      return 33;

    default:
        if (theKeySym <= 127){
            if (*pos < lvalue){
                if (*pos < 0) *pos = 0;     /* Restore default name */
                for (i=*pos; i<lvalue; i++) str[i]=' ';
                lsym = theKeySym;
                if (isalnum(lsym) || lsym == '/' || lsym == '.' || lsym == '_'){
                    str[*pos] = lsym;
                    *pos += 1;
                }
            }
        }
        break;
    }
    return 0;
}

/**********************************************************************/
int GetValue(XKeyEvent theKeyEvent, char str[], int lvalue, char tstri[], int *pos){
    int i;
    XComposeStatus theComposeStatus;
    KeySym theKeySym;
    int theBufferLength, theKeyBufferMaxLen = 4;
    char theKeyBuffer[5];

    theBufferLength = XLookupString(&theKeyEvent, theKeyBuffer,
        theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
    theKeyBuffer[theBufferLength] = '\0';

    switch(theKeySym){
    case XK_BackSpace:
    case XK_Delete:
          if (*pos > 0){
              for (i=*pos; i<lvalue; i++) str[i-1] = str[i];
              str[lvalue-1] = ' ';
          }
          *pos -= 1;
          break;

     case XK_question: return -2;
     case XK_Escape:   return -1;
     case XK_Return:
     case XK_KP_Enter: return 1;
     case XK_Tab:      return 2;
     case XK_R7:       return 11;
     case XK_Left:     return 21;
     case XK_Up:       return 12;
     case XK_Right:    return 23;
     case XK_Down:     return 32;
     case XK_R9:       return 13;
     case XK_R13:      return 31;
     case XK_R15:      return 33;

     default:
         if (*pos < lvalue){
             for (i=0; tstri[i] != 0; i++){
                 if (theKeySym == tstri[i]){
                     if (*pos < 0) *pos = 0;
                     str[*pos] = theKeyBuffer[0];
                     *pos += 1;
                 }
             }
         }
         break;
    }
    return 0;
}

/**********************************************************************/
int ufilebox_(INT_ *arr_size, char una[], char varNames[], char unad[]){
/* U-file name setting */
    Window theWindow;
    XEvent theEvent;
    const int wsym=8, hsym=13, /* symbol width and height */
        lvalue=10, lname=4,    /* length of value and name */
        xshif=5, yshif=3,      /* table corner and */
        lsym=4,
        n_columns=3;             /* # of columns */
    int UpLeftx, UpLefty,      /* window corner location */
        Width, Height,         /* window size */
        ix, iy, wbox, hbox,    /* parameter box corner and size*/
        nparam,                /* # of columns and parameters */
        ixa0=0, iya0=0, ixa, iya, /* arrow (textcursor) position */
        i, j, jbox, spos=0, vpos,
        iret, oldparam, ixold, iyold, ii=-1,
        icol, dcol, irow, drow;
    int lline, xButton, yButton, ind, ind1, ierr, key_flag;
    char Ufile_Name[40], vsym[5];

    strcpy(vsym, " => ");
    strcpy(Ufile_Name, "udb/");
    lline = 18;
    hbox = hsym+3;
    wbox = wsym*(lline+1);
    vpos = wsym*(lname+4);
    nparam = *arr_size;
    Width = 2*xshif + n_columns*wbox - wsym;
    Height =2*yshif + ((nparam - 1)/n_columns + 3)*hbox;
    jbox = 0;
    UpLeftx = 2;

    GetRWgeometry(&XRW, &YRW);
    UpLeftx = 2;
    j = XRW - Width - 8;
    UpLefty = YRW;
    if (j > UpLeftx) UpLeftx = j;
    j = YRW + Height + 30;
    if (j > theHeight) UpLefty = theHeight - Height - 30;
    theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0,
        "Save data in U-file format", 0,
        RootWindow(theDisplay, theScreen), theMenuCursor);
    XSelectInput(theDisplay, theWindow, POLL_EV_MASK);

    Change_Color(theGCA, 1, 0);  /* white background, black foreground */
    for (i=1; i <= nparam; i++){
        ix = wbox*((i - 1)%n_columns) + xshif;
        iy = hbox*((i - 1)/n_columns) + yshif;
        ind = (i - 1)*lname;
        WRITE_ theGCA, ix, hsym+iy, varNames+ind, lname);
        ind = (i - 1)*lvalue;
        for (j=0; j<lvalue; j++) una[ind+j] = unad[j];
        WRITE_ theGCA, ix+lname*wsym, hsym+iy, vsym, lsym);
        WRITE_ theGCA,ix+vpos, hsym+iy, una+ind, lvalue);
        XFlush(theDisplay);
    }
    oldparam = nparam;
    ixold = ix;
    iyold = iy;
    if (ii < 0) XWarpPointer(theDisplay, None, theWindow, 0, 0, 0, 0, xshif+40, yshif+13);
    Change_Color(theGCA, 50, 0);
    for (i=1; i < n_columns; i++){
        ind = xshif + i*wbox - wsym/2;
        LINGCA_ ind, 0, ind, Height - 2*hbox);
    }
    i = Height-2*hbox;
    LINGCA_ 0, i, Width, i);
    i--;
    LINGCA_ 0, i, Width, i);
    Change_Color(hghGC, 50, 0);
    WRITE_ hghGC, xshif, Height-4, " OK ", 4);
    XDrawRectangle(theDisplay, theWindow, hghGC, 5L, Height-17L, 30L, 16L);
    Change_Color(theGCA, 1, 0);
    i = n_columns*wbox/wsym-5;
    WRITE_ theGCA, xshif+4*wsym, Height-4, "/<ESC> - done;    Button, <TAB> or Arrow - select   ", i);
    WRITE_ theGCA, xshif, Height-4-hbox, "           Select box and enter U-file name              ", i+4);

// Table_control:
    while(1){
        XNextEvent(theDisplay, &theEvent);
        if (theEvent.xany.window == theRootWindow ){
            ProcessRootWindowEvent(&theEvent);
            continue;
        }
        key_flag = 0;
        if (theEvent.type == ButtonPress){
            xButton = theEvent.xbutton.x;
            yButton = theEvent.xbutton.y;
            if (jbox){
                oldparam = jbox;
                ixold = ix;
                iyold = iy;
            }
            jbox = FindBoxNum(xButton-xshif+wsym/2, yButton-yshif, wbox, hbox, n_columns, nparam);
            i = FindBoxNum(xButton-xshif, yButton-Height+hbox+1, 4*wsym, hbox, 1, 1);
            if (i){
                iret = -1;
                jbox = oldparam;
                ixold = ix;
                iyold = iy;
                key_flag = 1;
            }
            else{
                if (jbox == 0) continue;
                ind = (oldparam - 1)*lvalue;
                i = jbox;
                ix = wbox*((i - 1)%n_columns) + xshif;
                iy = hbox*((i - 1)/n_columns) + yshif;
                ind = (oldparam - 1)*lname;
                ind1 = (i - 1)*lname;
                WRITE_ theGCA, ixold, hsym+iyold, varNames+ind, lname);
                WRITE_ hghGC, ix, hsym+iy, varNames+ind1, lname);
                WRITE_ theGCA, ixold+lname*wsym, hsym+iyold, vsym, lsym);
                WRITE_ hghGC, ix+lname*wsym, hsym+iy, vsym, lsym);
                ind = (oldparam - 1)*lvalue;
                ind1 = (i - 1)*lvalue;
                WRITE_ theGCA, ixold+vpos, hsym+iyold, una+ind, lvalue);
                WRITE_ hghGC, ix+vpos, hsym+iy, una+ind1, lvalue);
                spos = 0;
                ixa = ix + vpos;
                iya = iy + hsym + 2;
                MoveArrow(theWindow, ixa0, iya0, ixa, iya);
                ixa0 = ixa;
                iya0 = iya;
                XFlush(theDisplay);
                continue;
            }
        } // end if "button press"
        else{
            if (jbox == 0){
                if (GetEsc(theEvent.xkey)){
                    break; // destroy u-file menu window
                }
                else{
                    continue;
                }
            }
            else{
                if (theEvent.type == KeyPress){
                    ind = (jbox-1)*lvalue;
                    iret = GetName(theEvent.xkey, una+ind, lvalue, &spos);
                    key_flag = 1;
                }
            }
        }

        if (key_flag == 1){
            if (spos >= 0){
                ixa = ix + vpos + spos*wsym;
                iya = iy + hsym + 2;
                MoveArrow(theWindow, ixa0, iya0, ixa, iya);
                ixa0 = ixa;
                iya0 = iya;
            }
            else{
                spos = -1;
                 strncpy(una+ind,unad,lvalue);
            }
            i = n_columns*wbox/wsym - 5;
            WRITE_ theGCA, xshif+4*wsym, Height-4,
                "/<ESC> - done;    Button, <TAB> or Arrow - select   ", i);
            if (strncmp(una+ind, unad, lvalue)){
                stcopy(Ufile_Name+4, una+ind, spos);
                ierr = 0;
                ierr = access(Ufile_Name, F_OK);
                if (!(ierr)){
                    if ( !(strncmp(una+ind, "          ", lvalue))){
                        strncpy(una+ind, unad, lvalue);
                    }
                    else{
                        WRITE_ hghGC, xshif+19*wsym, Height-4, "              ", 14);
                        WRITE_ hghGC, xshif+33*wsym, Height-4, " - file already exists", 22);
                        WRITE_ hghGC, xshif+(33-4-spos)*wsym, Height-4, Ufile_Name, 4+spos);
                    }
                }
            }

            WRITE_ hghGC, ix+vpos, hsym+iy, una+ind, lvalue);
            if (iret ==  0) continue;
            if (iret == -1) break;
            oldparam = jbox;
            if (iret == 2){
                jbox++;
                if (jbox > nparam) jbox = 1;
            }
            else{
                stcopy(Ufile_Name+4, una+ind, lvalue);
                dcol = iret%10 - 2;
                icol = (jbox - 1)%n_columns + dcol;
                if (icol < 0) icol = 0;
                if (icol >= n_columns) icol--;
                drow = iret/10 - 2;
                irow = (jbox - 1)/n_columns + drow;
                if (irow < 0) irow = 0;
                if (irow > (nparam-1)/n_columns) irow--;
                jbox = irow*n_columns + icol + 1;
                if (jbox == oldparam && icol == n_columns-1 && drow != -1) jbox++;
                if (jbox == oldparam && icol == 0       && drow !=  1) jbox--;
                if (jbox < 1) jbox = 1;
                if (jbox > nparam) jbox = nparam;
            }
            ixold = ix;
            iyold = iy;
            spos = 0;
            i = jbox;
            ix = wbox*((i - 1)%n_columns) + xshif;
            iy = hbox*((i - 1)/n_columns) + yshif;
            ind = (oldparam - 1)*lname;
            ind1 = (i - 1)*lname;
            WRITE_ theGCA, ixold, hsym+iyold, varNames+ind, lname);
            WRITE_ hghGC, ix, hsym+iy, varNames+ind1, lname);
            WRITE_ theGCA, ixold+lname*wsym, hsym+iyold, vsym, lsym);
            WRITE_ hghGC, ix+lname*wsym, hsym+iy, vsym, lsym);
            ind = (oldparam - 1)*lvalue;
            ind1 = (i - 1)*lvalue;
            WRITE_ theGCA, ixold+vpos, hsym+iyold, una+ind, lvalue);
            WRITE_ hghGC, ix+vpos, hsym+iy, una+ind1, lvalue);
            spos = 0;
            ixa = ix + vpos;
            iya = iy + hsym + 2;
            MoveArrow(theWindow, ixa0, iya0, ixa, iya);
            ixa0 = ixa;
            iya0 = iya;
            XFlush(theDisplay);
        } // if key_flag
    } // end Table_control;

    XWarpPointer(theDisplay, None, theRootWindow, 0, 0, 0, 0, Kevent.xcur, Kevent.ycur);
    XFlush(theDisplay);
    XDestroyWindow(theDisplay, theWindow);
    ibcursor = -1;
    XFlush(theDisplay);
    return 0;
}

/**********************************************************************/
int FindBoxNum(int ix, int iy, int wBox, int hBox, int n_columns, int nBox){
/* Returns box # or 0
  ix, iy:          current coordinates
  n_columns, nBox:     #columns, #boxes
  wBox, hBox:      box width, height
*/
    int icol, pnum;

    if (ix <= 0 || iy < 0) return 0;
    icol = ix/wBox;
    if (icol >= n_columns) return 0;
    pnum = (iy/hBox)*n_columns + icol + 1;
    if (pnum > nBox) return 0;
    return pnum;
}
