#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
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
#define RETURNPOINTER_ XWarpPointer(theDisplay, None, theRootWindow, 0, 0, 0, 0, Kevent.xcur, Kevent.ycur)
#define LINGCA_ XDrawLine(theDisplay, theWindow, theGCA,

#define F1sw 8      /*font 1 symbol width */
#define F1sh 13     /*font 1 symbol height */
#define BTyof 0     /*button text vertical offset */
#define BTxof 0     /*button text horizontal offset */
#define BVers 4     /*button vertical separation */
#define NMRMB 36    /*main Review menu table */
int nmamb = NMRMB;  /*main Astra  menu table */
char *MAMT[NMRMB];  /* Dimension is defined as */
char MAMK[NMRMB];   /* max(nmamb,NMRMB) */
#define Button struct BUTTON
Button {int mexl, meyd, mexr, meyu, melen; char keysym, *mename;};
Button MAMB[NMRMB];
#define NDW 12      /* number of dialog windows */
char *DWTitle[NDW] = {
    "Variable control", /*No. 1  "V"*/
    "Constant control", /*No. 2  "C"*/
    "Times & Grids",    /*No. 3  "D"*/
    "Window control",   /*No. 4  "W"*/
    "Scale control",    /*No. 5  "S"*/
    "Y-shift",          /*No. 6  "Y"*/
    "Time interval",    /*No. 7  "M" in the mode 7*/
    "Mark times:  < 0 - skip,  0 - dim,  > 0 - color #", /*"M" modes 4,5*/
    "Equilibrium control",   /*No. 9  call from metric */
    "1D_Ufile",              /*No. 10 Not used */
    "2D_Ufile",              /*No. 11 Not used */
    "NBI const for beam No", /*No. 12 call from nbiext */
};
int DWnamlen[NDW] = {
    6, /*No. 1  PRNAME*/
    6, /*No. 2  CFNAME*/
    6, /*No. 3  DTNAME*/
    4, /*No. 4  NAME[TR]*/
    4, /*No. 5  NAME[TR]*/
    4, /*No. 6  NAME[TR]*/
    6, /*No. 7  NAM7*/
    6, /*No. 8  NAMEP*/
    6, /*No. 9  DTNAME*/
    6, /*No. 10 1D_Ufi*/
    6, /*No. 11 2D_Ufi*/
    6, /*No. 12 NBI control*/
};

struct {int xcur, ycur;} Kevent;

void Change_Color(GC, int, int);
void Menu_table(int, int, int, char[], char*[], Button[], int);
void Put_button(Button, GC);
void num2str(double, char*, int);
void MoveArrow(Window, int, int, int, int);
void changeGCcolor(GC, INT_*);
void PutColorName(Window, int, int, int, int);
void stcopy(char*, char*, int);
int isascii(int);
int isprint(int);
int isalnum(int);
int nextevent(INT_*, INT_*, INT_*, int, Button[], char[], char*[]);
int asklis_(INT_*, double*, char[], INT_*);
int askcol_(char[], char[], char[], INT_*, INT_*, INT_*, INT_*);
int asktab_(char[], char[], char[], INT_*, INT_*, INT_*, INT_*);
int askgrf_(char[], char[], char[], INT_*, INT_*, INT_*, INT_*, INT_*);
int ufilebox_(INT_*, char[], char[], char[]);
int FindBoxNum(int, int, int, int, int, int);
int GetEsc(XKeyEvent);
int GetValue(XKeyEvent, char[], int, char[], int*);
int warning_dialogbox_();
int If_empty(int, int, int, int[], int[], char[], int);
int GetKey(XKeyEvent, char[], int*);
Window Open_Window(int, int, int, int, int, char[], int, Window, Cursor);

/**********************************************************************/
void xaxis(INT_ *modex){
    xaxis_(modex);
}

void xaxis_(INT_ *modex){ /* Change of X-coordinate for plots */
    switch (*modex){
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
    int Xx, Xy, i, ixx=8, dx=5, iyy;
    int BHors=5;                /*button horizontal separation */
    char *BUTEXT[NMRMB] =
        {"16*f(a)", "8*f(a)", "Refresh", "2*f(a,t)",
         "8*f(t)", "User graph", "2*f(R,t)", "Phase space",
         "Next", "8*f(psi)", "Equil", "Backward",
         "Scales", "Variables", "Type data", "Port_PS",
         "Windows", "Constants", "Save log" , "Land_PS",
         "Select", "Grids", "Write data", "U-files",
         "Style", "Type model", "What X-axis", "Y-shift",
         "Run", "Step", "Quit", "Help"};
    char BUTKEY[NMRMB] = {
        '1', '2', 'R', '4', '6', '9', '5', '7', 'N', '3', '8', 'B',
        'S', 'V', 'T', 'G', 'W', 'C', 'I', 'Q', 'M', 'D', 'F', 'U',
        '.', 'L', 'X', 'Y', '\015', '\040', '\057', 'H'};
    nmamb = 32;
/* Draw separating lines between plots and menu */
    Change_Color(theGCA, AstraColorNum[2], AstraColorNum[3]);
    XDrawLine(theDisplay, theRootWindow, theGCA, 0, XWH-128, XWW-1, XWH-128);
    XDrawLine(theDisplay, theRootWindow, theGCA, 0, XWH-110, XWW-1, XWH-110);
    XDrawLine(theDisplay, theRootWindow, theGCA, 0, XWH-109, XWW-1, XWH-109);
    for (i=0; i<nmamb; i++){
        MAMT[i] = BUTEXT[i];
        MAMK[i] = BUTKEY[i];
    }
    xaxis(modex);
    iyy = XWH - 74;
    for (i=0; i<8; i++){
        Menu_table(ixx + i*(XWW/8+dx), iyy, 4, MAMK+4*i, MAMT+4*i, MAMB+4*i, BHors);
    }
    if (ibcursor >= 0) Put_button(MAMB[ibcursor], hgh_menuGC);
/* Menu titles */
    Change_Color(hghGC, AstraColorNum[2], AstraColorNum[3]);
    Xx = ixx + 40;
    Xy = iyy - 5;
    XDrawImageString(theDisplay, theRootWindow, hghGC, Xx, Xy,
        "Graphic mode           Presentation       Control", 49);
    Xx += 475;
    XDrawImageString(theDisplay, theRootWindow, hghGC, Xx, Xy, "In/Out     Status", 17);
}

/**********************************************************************/
void Menu_table(int xm, int ym, int nbutt, char mek[], char *met[], Button butt[], int BHors){
    int i, leng, lenm, xwid;

    lenm = 0;
    for (i=0; i<nbutt; i++){
        leng = strlen(met[i]);
        if (leng > lenm) lenm = leng;
    }
    xwid = lenm*F1sw + 2*BTxof + BHors;
    for (i=0; i<nbutt; i++){
        leng = strlen(met[i]);
        butt[i].mexl = xm;
        butt[i].mexr = butt[i].mexl + xwid;
        butt[i].meyu = ym + i*(F1sh + 2*BTyof + 1+BVers);
        butt[i].meyd = butt[i].meyu + F1sh + 2*BTyof + 1;
        butt[i].melen  = leng;
        butt[i].mename = met[i];
        butt[i].keysym = mek[i];
    }
    for (i=0; i<nbutt; i++){
        Put_button(butt[i], theGCA);
    }
}

/**********************************************************************/
void Put_button(Button but, GC aGC){
/* Button drawing */
    int yd, xd, yu, xu;
    yd = but.meyd;
    xd = but.mexl;
    yu = but.meyu;
    xu = but.mexr;
    Change_Color(aGC, 1, 0);
    XClearArea(theDisplay, theRootWindow, xd-2, yu-2, xu-xd+4, yd-yu+4, False);
    XDrawImageString(theDisplay,theRootWindow, aGC, xd+BTxof+2, yd-BTyof-3, but.mename, but.melen);

    XDrawLine(theDisplay, theRootWindow, aGC, xd+2, yd  , xu-2, yd  );
    XDrawLine(theDisplay, theRootWindow, aGC, xu-2, yd  , xu  , yd-2);
    XDrawLine(theDisplay, theRootWindow, aGC, xu  , yd-2, xu  , yu+2);
    XDrawLine(theDisplay, theRootWindow, aGC, xu  , yu+2, xu-2, yu  );
    XDrawLine(theDisplay, theRootWindow, aGC, xu-2, yu  , xd+2, yu  );
    XDrawLine(theDisplay, theRootWindow, aGC, xd+2, yu  , xd  , yu+2);
    XDrawLine(theDisplay, theRootWindow, aGC, xd  , yu+2, xd  , yd-2);
    XDrawLine(theDisplay, theRootWindow, aGC, xd  , yd-2, xd+2, yd  );
/* No round caps:
    XDrawLine(theDisplay, theRootWindow, aGC, xd, yd, xu, yd);
    XDrawLine(theDisplay, theRootWindow, aGC, xu, yd, xu, yu);
    XDrawLine(theDisplay, theRootWindow, aGC, xu, yu, xd, yu);
    XDrawLine(theDisplay, theRootWindow, aGC, xd, yu, xd, yd);
*/
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
        if (i >= 0 && i <= nmamb) *key = (int)MAMK[i];
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
    i = nextevent(theKey, xCursor, yCursor, nmamb, MAMB, MAMK, MAMT);
    while (i == 0 && *theKey == 0 ){
        i = nextevent(theKey, xCursor, yCursor, nmamb, MAMB, MAMK, MAMT);
    }
    return i;
}

/*********************** Waiting events *******************************/
int nextevent(INT_ *theKey, INT_ *xCursor, INT_ *yCursor, int nmbt,
    Button mbb[], char mbk[], char *mbt[]){
/* ibcursor -     current box # or -1
   iact     - highlighted box # or -1 */
    XEvent theEvent;
    XKeyEvent theKeyEvent;
    XComposeStatus theComposeStatus;
    KeySym theKeySym;
    int length, ibox, longKey, theKeyBufferMaxLen=4;
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
        ibox = Cursor_in_Box();
        if (ibox >= 0) *theKey = (int)mbk[ibox];
        if (theDepth == 1 && *theKey == 'A') *theKey = 0;        /*No color table*/
        if (*theKey == 015) Put_button(mbb[ibox], theGCA);  /*Undo highlight*/
        if (*theKey != 0) ibox = -1;
        switch(ibox){
        case(-1): return 65001;
        case(24): return 65363;  /*   >    */
        case(25): return 65502;  /*  End   */
        case(26): return 65361;  /*   <    */
        case(27): return 65496;  /*  Home  */
        case(28): return 65364;  /*   >>   */
        case(29): return 65504;  /*  PgDn  */
        case(30): return 65362;  /*   <<   */
        case(31): return 65498;  /*  PgUp  */
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
    for (i=0; i<nmamb; i++){
        if (Kevent.xcur >= MAMB[i].mexl && Kevent.xcur <= MAMB[i].mexr &&
            Kevent.ycur >= MAMB[i].meyu && Kevent.ycur <= MAMB[i].meyd){
            return i;
        }
    }
    return -1;
}

/**********************************************************************/
void mvcursor_(INT_ *key, INT_ *ix, INT_ *iy){
/* Move cursor by one pixel */
    switch(*key){
    case(361):
        *ix--; /* Left  */
        break;
    case(362):
        *iy--; /* Up    */
        break;
    case(363):
        *ix++; /* Right */
        break;
    case(364):
        *iy++; /* Down  */
        break;
    }
    XWarpPointer(theDisplay, None, theRootWindow, 0, 0, 0, 0, *ix, *iy);
}

/**********************************************************************/
int asklis_(nofbox, array, theNames, id)
    INT_ *nofbox, *id;
    char theNames[];
    double *array;
{        Window                theWindow;
        XEvent                theEvent;
        double                valn;
        int        UpLeftx =2, UpLefty =10,/* window corner location        */
                Width =452, Height =480;/* window size                        */
        int        ix, iy, wbox, hbox,         /* parameter box corner and size*/
                wsym = 8, hsym = 13,        /* symbol width and height        */
                lvalue = 6, lname = 6,        /* length of value and name        */
                xshif = 5, yshif = 3,        /* table corner                 */
                nclmn = 4, nparam,         /* # of columns and parameters        */
                spos  = 0, vpos,         /* edit & edit start positions        */
                ixa0=0, iya0=0, ixa,iya,/* arrow (textcursor) position        */
                namlen, ibox, iret, oldparam, ixold, iyold;
        int        lline, xButton, yButton, i, ii=-1, ind, ind1,
                ihelp, selalb=0,  contr, icol, dcol, irow, drow;
        float        param;
        char        value[10],  ovalue[10],  stri[10],  vsym = '=',title[70],
                exp_name[60], equ_name[60],
                stri50[50],stri40[40], stri256[256], theName[10], legend[128];
        i = *id-1;        namlen = DWnamlen[i];        strcpy(title,*(DWTitle+i));

        lline = lname + lvalue + 1;
        if ((lvalue == 0) || (lname == 0)) lline = lname+lvalue;
        hbox = hsym+3;                wbox = wsym*(lline+1);        vpos= wsym*(lname+1);
         selalb = 1;        nparam   = *nofbox;        contr=1;
        if( nparam==0 ) return (-1);
        if( nparam<0 ) { nparam = -nparam; selalb = 0; }
        if( nparam>1000 ) { nparam=nparam-1000;        contr=0; }
        if( *id == 3 )  contr=1;
        Width =2*xshif+nclmn*wbox-wsym;
        Height =2*yshif+((nparam-1)/nclmn+2)*hbox;
        ibox = 1;        if (contr == 0 ) ibox = 0;
               if( *id == 8 )  { ibox = 0;        selalb = -1; }

        GetRWgeometry (&XRW,&YRW);
        UpLeftx = 2;     i = XRW-Width-8;                UpLefty = YRW;
        if (i > UpLeftx)  UpLeftx = i;                i = YRW+Height+30;
        if (i > theHeight) UpLefty = theHeight-Height-30;
        theWindow
          = Open_Window (UpLeftx, UpLefty, Width, Height, 0, title, 0,
             RootWindow(theDisplay,theScreen), theMenuCursor);
        XSelectInput (theDisplay, theWindow, POLL_EV_MASK);        ihelp = 0;
Create_table:
        Change_Color (theGCA, 1, 0);        /* white background, black foreground */
        for (i=1; i <= nparam; i++)         {
        ix =wbox*((i-1)%nclmn)+xshif;        iy =hbox*((i-1)/nclmn)+yshif;
        if ( lname ) {        ind = (i-1)*namlen;
                WRITE_ theGCA, ix, hsym+iy, theNames+ind, namlen); }
        if ( lname*lvalue ) {
                WRITE_ theGCA, ix+vpos-wsym, hsym+iy, &vsym, 1); }
        if ( lvalue ) { valn=*(array+i-1); num2str(valn, value, lvalue);
                WRITE_ theGCA, ix+vpos, hsym+iy, value, lvalue); }
        XFlush(theDisplay);                }
        oldparam = nparam;        ixold = ix;        iyold = iy;
        for (i=0; i<lvalue; i++)         ovalue[i]=value[i];
        if (ibox == 1 && selalb != -1 && ii < 0) MVPOINTER_ xshif+40,yshif+13);
        Change_Color (theGCA, 50, 0);
        for (i=1; i < nclmn; i++) { ind = xshif+i*wbox-wsym/2;
        LINGCA_ ind,0,ind,Height-hbox); }
        i = Height-hbox; LINGCA_ 0,i,Width,i);
        i--;                  LINGCA_ 0,i,Width,i);
        WRITE_ hghGC, xshif, Height-3, "OK",2);
        if ( ibox == 0 && contr == 0 ) MVPOINTER_ xshif+5,Height-5);
               if ( *id == 8 )  MVPOINTER_ xshif+Width-50,Height-5);
        i = nclmn*wbox/wsym-5;        ind = 50;        if(i<50)        ind=i;
/*      Modified by Antoine Merle on 09/18/2013
        I have modified the layout of the legend bar to make it fit into
        the 70 character length corresponding to its initialization, since
        this could cause a code crash due to buffer overflow on some systems */
        strcpy(legend,"/ Esc - done;");
        if (*id <= 6) { strcat(legend," Button/Tab - select;");
                        strcat(legend," Return/Tab - enter;");
                        strcat(legend," ? - quick help"); }
        else          { strcat(legend,"     Button/Tab - select;");
                        strcat(legend,"     Return/Tab - enter;"); }
        Change_Color (hintGC, 42, 0);        ind = strlen(legend);
        if (contr == 1) WRITE_ hintGC, xshif+2*wsym, Height-3, legend, ind);
        Change_Color (hintGC, 1, 0);
        if (contr == 0) WRITE_ theGCA, xshif+4*wsym, Height-3,
        "    Information table.   No changes permitted.     ", ind);
        if ( *id == 8 ) {
              if (selalb == 0) WRITE_ hintGC,
                        xshif+8+(i-6)*wsym,Height-3," (Select all ) ",15);
              if (selalb == -1) WRITE_ hintGC,
                        xshif+8+(i-6)*wsym,Height-3,"(Unselect 0)",12);
                        }
        Change_Color (theGCA, 1, 0);
        if (ibox) goto Newparam;
Table_control:
        XNextEvent (theDisplay, &theEvent);
        if ( theEvent.xany.window == theRootWindow )
           { /*  printf("theRootWindowEvent.type =%d\n",theEvent.type); */
             ProcessRootWindowEvent (&theEvent);
             goto Table_control;
           }
        if(theEvent.type == Expose) { ii = 0;  goto Create_table;}
        if(theEvent.type == ButtonPress)
           {	xButton	= theEvent.xbutton.x;
		yButton = theEvent.xbutton.y;
		ind  = FindBoxNum(xButton-xshif, yButton-Height+hbox+1,
						 	4*wsym, hbox, 1, 1);
		if(contr)
		{
		if (ibox) { oldparam = ibox;	ixold = ix;	iyold = iy;
			if (spos>0)  { 	sscanf(stri,"%6g",&param);
					*(array+ibox-1)=param; valn=param;
					num2str(valn, value, lvalue); spos=-1;
				     }
			for (i=0; i<lvalue; i++) 	ovalue[i]=value[i];
			  }
		ibox = FindBoxNum(xButton-xshif+wsym/2, yButton-yshif,
				  wbox, hbox, nclmn, nparam);
		ind1 = FindBoxNum(xButton-xshif-nclmn*wbox+10*wsym,
				yButton-Height+hbox+1, 10*wsym, hbox, 1, 1);
		if(ind1!=0 && selalb!=1)
		   { for (i=0; i<nparam; i++)
			{ if( *(array+i)<=0 ) *(array+i)=selalb; }
			{ if(selalb)  selalb=0; else selalb=-1; }
			  goto Create_table;
		   }
		}
		if( ind ) { iret = -1;	ibox = oldparam;  goto Escend; }
		if (ibox == 0)	goto Table_control;
Newparam:
		ix =wbox*((ibox-1)%nclmn)+xshif;
		iy =hbox*((ibox-1)/nclmn)+yshif;
		if ( lname ) {ind = (oldparam-1)*namlen; ind1 = (ibox-1)*namlen;
		    WRITE_ theGCA, ixold, hsym+iyold, theNames+ind, namlen);
		    Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
 		    WRITE_ hghGC, ix, hsym+iy, theNames+ind1, namlen);
		    Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]); }
		if ( lname*lvalue ) {
		    WRITE_ theGCA, ixold+vpos-wsym, hsym+iyold, &vsym, 1);
		    Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
		    WRITE_ hghGC, ix+vpos-wsym, hsym+iy, &vsym, 1);
		    Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]); }
		if ( lvalue ) {valn=*(array+ibox-1); num2str(valn,value,lvalue);
		    WRITE_ theGCA, ixold+vpos, hsym+iyold, ovalue, lvalue);
		    Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
		    WRITE_ hghGC, ix+vpos, hsym+iy, value, lvalue);
		    Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]);
		    spos = -1; for (i=0; i<lvalue; i++) stri[i]=' ';
		    ixa = ix+vpos;	iya = iy+hsym+2;
		    MoveArrow(theWindow,ixa0,iya0,ixa,iya);
		    ixa0 = ixa;		iya0 = iya; }
		XFlush(theDisplay);
		goto Table_control;
	   }
	if ((ibox == 0)||(lvalue == 0)) {
			if(GetEsc(theEvent.xkey)) goto EndDialog;
				else	goto Table_control; }
	iret = GetValue(theEvent.xkey, stri, lvalue, "0123456789.-+eE", &spos);
Escend:	if ( contr==0 ){/*printf("iret %d\n",iret);*/	goto EndDialog;}
	if ( spos >= 0 )
	   {	ixa = ix+vpos+spos*wsym;	iya = iy+hsym+2;
		MoveArrow(theWindow,ixa0,iya0,ixa,iya);
	        ixa0 = ixa;			iya0 = iya; }
	if ( iret && (spos>0) )
	   { 	sscanf(stri,"%6g",&param);
		/*printf("input =%s,    double =%g \n",stri,param);*/
		*(array+ibox-1) = param; valn=param;
		num2str(valn, value, lvalue); spos=-1;	}
	if ( iret == -1 )  	goto EndDialog;
	if ( iret == -2 )
	   { i = 0;	ii = namlen;
TableQuestionMark:    strncpy(theName,theNames+(ibox-1)*namlen+i,ii);
	     if ( theName[0]  == ' ' )	{ii--;	i++;  goto TableQuestionMark;}
	     theName[ii] = '\0';
	     while ( theName[--ii] == ' ' ) theName[ii] = '\0';	ii++;

	     if (*id == 1)		/*Variables <- for/const.in*/
		{ strcpy(stri50,"grep -i \"");   strncat(stri50,theName,ii);
	          strcat(stri50," \" main/variables.txt");   system(stri50);  }
	     if (*id == 2)		/*Constant table*/
		{
		  strcpy (stri256,"grep -i -w ");
		  strncat(stri256,theName,ii);
   		  strcat (stri256," tmp/model.tmp");
		  ii = system(stri256);
		  if (ii == 256) {
                      printf("The constant %s is not assigned\n", theName);
		  }
		  }
	     if (*id == 3)		/*Time, grid control*/
		{  strcpy(stri50,"grep -i \"");  strncat(stri50,theName,ii);
	           strcat(stri50," \" main/internal.txt");   system(stri50);  }
	     if (*id == 4)
		{  strcpy (stri256,"grep -w ");    strncat(stri256,theName,ii);
		   strcat (stri256," tmp/model.txt");
		   if (!ihelp) printf(" No. Scale  Name  Output expression\n");
		   system(stri256);	ihelp = 1;	}
	     if (*id == 5)		/*Scale control*/
		{  strcpy (stri256,"grep -w ");    strncat(stri256,theName,ii);
		   strcat (stri256," tmp/model.txt");
		   if (!ihelp) printf(" No. Scale  Name  Output expression\n");
		   system(stri256);	ihelp = 1;	}
	     if (*id == 6)		/*Vertical shift*/
		{  strcpy (stri256,"grep -w ");    strncat(stri256,theName,ii);
		   strcat (stri256," tmp/model.txt");
		   if (!ihelp) printf(" No. Scale  Name  Output expression\n");
		   system(stri256);	ihelp = 1;	}
	     if (*id == 7)
		{  printf("Curve presentation\n");	}
	     if (*id == 8)
		{  printf("Mark time slices\n");	}
		/* Control print:	printf("[%s]\n",stri50);  */
	   }
	if ( spos < 0 )
	   {	spos = -1; for (i=0; i<lvalue; i++) stri[i]=' ';
		Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
		WRITE_ hghGC, ix+vpos, hsym+iy, value, lvalue);
		Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]); }
	else
	   {    Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
		WRITE_ hghGC, ix+vpos, hsym+iy, stri, lvalue);
		Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]); }
	if ( iret == 2 )
	   {	oldparam = ibox;	ixold = ix;	iyold = iy;
		for (i=0; i<lvalue; i++) 	ovalue[i]=value[i];
		ibox++; if(ibox > nparam) ibox=1; goto Newparam; }
	if ( iret > 2 )
	   {	oldparam = ibox;	ixold = ix;	iyold = iy;
		for (i=0; i<lvalue; i++) 	ovalue[i]=value[i];
 		dcol=iret%10-2; 	icol =(ibox-1)%nclmn+dcol;
		if (icol<0) icol=0;	if (icol>=nclmn) icol--;
		drow=iret/10-2;		irow =(ibox-1)/nclmn+drow;
		if (irow<0) irow=0;	if (irow>(nparam-1)/nclmn) irow--;
		ibox =irow*nclmn+icol+1;
		if(ibox == oldparam && icol == nclmn-1 && drow != -1) ibox++;
		if(ibox == oldparam && icol == 0       && drow != 1)  ibox--;
		if (ibox < 1) ibox=1;	if (ibox > nparam) ibox=nparam;
		goto Newparam; }
				goto Table_control;
EndDialog:
	RETURNPOINTER_;
	XSetInputFocus(theDisplay,theRootWindow,
	         RevertToPointerRoot,theEvent.xkey.time);
	XFlush(theDisplay);
	XDestroyWindow(theDisplay, theWindow);
	ibcursor = -1;
	XFlush(theDisplay);
	return (0);
}

/**********************************************************************/
int asktab_(title, template, array, len, nrows, ngroup, morow)
     INT_     *len, *nrows, *ngroup, *morow;
     char    title[], template[], array[];
/*
Called from ASXWIN and ASTWIN (file surv.f90) that in turn are invoked by key "M" from IFKEY

Input:	title	- Title of the table
	template  string defining a structure of the table and its 1st line
		  1st line does not appear if all non-'|' symbols are spaces
	array	- data (numbers or strings) for input and output
	len	- length of "array" element according to description
			in the calling routine (80 in the example below)
	nrows	- number of rows in a table to be created
	ngroup	- if > 0 distance (in rows) between blue separating lines
	morow	- if > 0 separates bottom of the table with a fat blue line

Returned value:
	-1  - Error (window was not created)
	 0  - Normal exit (table updated)
	>0  - Temporary exit (the returned value goes to the calling routine)
Call from FORTRAN
    Example 1:
	integer		nrow, len
	character*80	title, template, array(10)
	len = 80
	nrow = 5	! must be <= 10
	title = "Example"//char(0)
	template = " field1 |field2|    | f4 |       "//char(0)
	array(1) = " CumBol &  10.1& 2  &    &-1.E+5 "
	array(2) = "                                 "	! Empty string
	array(3) = "                                 "	! Empty string
	array(4) = "                                 "	! Empty string
	array(5) = "                                 "	! Empty string
C symbols in the positions "|" will not appear in the table
C 4th parameter must be 80
C 5th parameter must be <= 10
	call	asktab(title,template,array,80,3,0,0)

More examples in the file "for/surv.f" subroutine "ASXWIN" or "ASTWIN"

Add features:
   If (Width=2*xshif+wsym*(last_pos) < 37*wsym)
       or (Height=2*yshif+(nlines+3)*hbox > theHeight) double Width
*/
{
  XEvent	theEvent;
  static	Window	theWindow;
  static	int	xshif = 5, yshif = 3,	/* table corner	*/
    ind, ibold, hbox,	 /* current/old box No. & height	*/
    irow0=0, 	wbox, 	      /* number of chars in box 	*/
    Width =452, Height =480,     /* defaut window sizes		*/
    inold, ixold, iyold, icold, ixa0, iya0, arrdim,
    nclmn, nwid[20], nsta[20], mode=0, ii=0,
    icol, isym, irow, ixa, iya,  /* arrow (textcursor) position */
    spos=0,	   /* abs. and rel. position of symbol in table */
    lvalue,		/* 1 if box was changed, 0 otherwise	*/
    wsym = 8, hsym = 13;	/* symbol width and height	*/
  int	UpLeftx, UpLefty =10,	/* window corner location	*/
    nend, ibox, iret, ix, ix1, iy, i, j,
    mxfields = 20;		/* max / actual # of columns	*/
                  /* int      count, mode = QueuedAfterReading; */
  static	char	stri[132];
  if ( mode == 1 )  goto	Table_control;		mode = 0;
  j = strlen(title);	if ( j > 132 )			{
    printf("%s\n\"%s\"\n%s\n","ASKTAB >>> Title",title,
	   "           is too long");		return -1;	}
  if  ( *len > 132 )				{
    printf("%s \"%s\" %s\n","ASKTAB >>> Table",title,
	   ". Input string is too long");	return -1;	}
	nend = strlen(template);
	if ( nend > 132 )	{
		printf("%s \"%s\" %s\n","ASKTAB >>> Table",title,
		     ". Requested width is too large");	return -1;	}
	for ( j=0; j < nend ; j++)
	      if ( template[j] != '|'  &&  template[j] != ' ' ) irow0 = 1;

		/****   Split "template" in blocks  ****/
	for ( j=nclmn=nsta[0]=0; j < nend ; j++ )
	    { if ( template[j] == '|' )
		 { nwid[nclmn] = j-nsta[nclmn]; nsta[++nclmn] = j+1;}
	      if ( nclmn > mxfields ) {
		printf("%s %d\n","ASKTAB >>>  Too many input fields ",nclmn);
			return -1;	}
	    }
	nwid[nclmn] = nend-nsta[nclmn];		nclmn++;	arrdim = *len;

	ixa0  = ibold = 0;	UpLeftx = 2;   	hbox   = hsym+3;
	Width = 2*xshif+nend*wsym;   Height = 2*yshif+(*nrows+4+irow0)*hbox;
	if (Kevent.xcur+Kevent.ycur < 10)
	  { Kevent.xcur = 330;
	    Kevent.ycur = 430;
	  }
       	/*  j = XWX-Width-10;	if (j > UpLeftx)  UpLeftx = j; */
	GetRWgeometry (&XRW,&YRW);
	UpLeftx = 2;	   j = XRW-Width-8;	UpLefty = YRW;
	if (j > UpLeftx)   UpLeftx = j;		j = YRW+Height+30;
	if (j > theHeight) UpLefty = theHeight-Height-30;
	theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0, title, 0,
		    RootWindow(theDisplay,theScreen), theMenuCursor);
	MVPOINTER_ xshif+wsym*nwid[0]/2,yshif+irow0*hbox+hsym);
	XSelectInput (theDisplay, theWindow, POLL_EV_MASK);

Create_table:
	Change_Color (theGCA, 3, 0); /* white background, blue  foreground */
	if  ( irow0 )
	    {	/*********  Type template in blue  *********/
	    for ( j=0, ix=xshif; j < nclmn ; j++)
		{ WRITE_ theGCA, ix, hsym, template+nsta[j], nwid[j]);
		  if (j+1 == nclmn) break;	ix += wsym*(nwid[j]+1);
		}
		/*********  Draw blue horizontal lines  *********/
	    j = hbox;			LINGCA_ 0,j+4,Width,j+4);
	    }
	if  ( *ngroup > 0 )
	    for ( j = (*ngroup+irow0)*hbox; j+4+yshif < Height-2*hbox ;
		  j += (*ngroup)*hbox )	  LINGCA_ 0,j+4,Width,j+4);
	if  ( *morow > 0 )
	    {	j = Height-(*morow+2)*hbox-3;	LINGCA_ 0,j,Width,j);
		++j;				LINGCA_ 0,j,Width,j);
	    }
		/*********  Draw blue vertical lines  *********/
	for ( j = 0, ix=xshif-0.5*wsym; j < nclmn-1; j++)
	    { ix += wsym*(nwid[j]+1);	ix1 = ix;
	      if (nwid[j] == 0)		ix1 = ix1-2;
	      if (nwid[j+1] == 0)	ix1 = ix1+2;
	      LINGCA_ ix1,0,ix1,Height-2*hbox);
	    }
		/*********  Draw bottom line & comments  *********/
	WRITE_ theGCA, xshif, Height-3-hbox,
		"Button - select,  <Tab>, <Ret>, Arrows - move", 45);
	WRITE_ theGCA, xshif+4*wsym, Height-4, "/<ESC> - done", 13);
	Change_Color (theGCA, 50, 0);  /* white background, red foreground */
	j = Height-2*hbox; 	LINGCA_ 0,j,Width,j);
	j++; 		   	LINGCA_ 0,j,Width,j);
	Change_Color (hghGC, 50, 0);
	WRITE_ hghGC,  xshif, Height-4, " OK ",4);
	XDrawRectangle (theDisplay,theWindow,hghGC,5L,Height-17L,30L,16L);
	Change_Color (theGCA, 1, 0);	Change_Color (hghGC, 1, 0);
					/* white background, black foreground */
		/*********  Draw input data  *********/
	for (   i=0; i < *nrows; i++ )				/* row loop */
 	    {	iy =yshif+hsym+hbox*(i+irow0);
	    for ( j=isym=0, ix=xshif; j < nclmn ; j++)	     /* column loop */
	    	{ if (j > 0) { isym = nsta[j]; ix += wsym*(nwid[j-1]+1); }
	    	  WRITE_ theGCA, ix, iy, array+i*arrdim+isym, nwid[j]);
		}
	    }
	if ( ibold > 0 )				/* after Expose event */
	   {	isym = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=xshif+wsym*isym;	iy=yshif+hsym+hbox*(irow+irow0-1);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	if ( ii == 0 )		/* Select Upper Left box on entry */
	   { ibox = ibold = irow = icol = icold =ii = 1;
		inold = ind = isym = lvalue = 0;	wbox = nwid[0];
		ix = ixold = ixa = ixa0 = xshif+wsym*isym;
		iy = iyold = iya = iya0 = yshif+hsym+hbox*irow0;
		WRITE_ hghGC, ix, iy, array, wbox);
		strncpy(stri,array,wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	XFlush(theDisplay);

Table_control:
	XNextEvent (theDisplay, &theEvent);
	if ( theEvent.xany.window == theRootWindow )
           { ProcessRootWindowEvent (&theEvent);
	     goto Table_control;
	   }
	if (theEvent.type == Expose)	goto	Create_table;
	if (theEvent.type == ButtonPress)
	{   ix = theEvent.xbutton.x;	iy = theEvent.xbutton.y;
	    irow = iy-yshif-hbox*irow0-2;
	    if ( irow >= 0)
	    {	irow = irow/hbox+1;		i=(ix-xshif)/wsym;
		if ( ibold > 0 )	strncpy(array+ind,stri,wbox);
		if ( irow > *nrows )		  /* No selection or Exit */
		   { if ( i > 3 || iy < Height-4-hsym ) goto Unselect;
			goto EndDialog;			/* OK was pressed */
		   }
		for ( j=0; j < nclmn; )			/* Selection is made */
		    { i -= nwid[j]+1;	icol=++j;   if ( i < 0 ) break;	}
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 )	goto	Unselect;
Fillbox:	isym = spos = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=ixa=xshif+wsym*isym;	iy=iya=yshif+hsym+hbox*(irow+irow0-1);
		ind = (irow-1)*arrdim+isym;	wbox = nwid[icol-1];
		ibox = (irow-1)*nclmn+icol;
		if ( ibold != 0 )
		     WRITE_ theGCA, ixold, iyold, stri, nwid[icold-1]);
		strncpy(stri,array+ind,wbox);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, ixa0, iya0, ixa, iya);
		icold = icol;	inold = ind;	ixold = ix;	iyold = iy;
		ixa0 = ixa;	iya0 = iya;	ibold = ibox;	lvalue = 0;
	    }
	    else
Unselect:   {	if ( ibold > 0 )			{
		WRITE_ theGCA, ixold, iyold, array+inold, nwid[icold-1]);
		MoveArrow(theWindow,ixa0,iya0,0,iy);	ibold=spos=0;	}
	    }
	    XFlush(theDisplay);
	    goto Table_control;
	}

	if ( ibold > 0 )
	{
	   iret = GetKey (theEvent.xkey, stri, &spos);
	   if ( iret == -1 && ibold > 0 )			/* <Esc> */
	      { strncpy(array+ind,stri,wbox);   goto EndDialog;
	      }
	   if ( iret == 0 )
	      {	if (spos == 0 ) { 		/* 1-st entry in the box */
		for (j=1; j < 132; j++) stri[j]=' '; stri[wbox] = '\0';  }
		if (++spos > wbox) spos = wbox;		lvalue = 1;
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 1 )				/* <Del> or <BS> */
	      {	if (--spos >= 0) { strncpy(stri+spos,stri+spos+1,wbox-spos);
				  stri[wbox-1] = ' ';	}
		if (spos < 0) { spos = 0;  strncpy(stri,array+ind,wbox); }
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( (iret == 2 || iret == 3) && ibold > 0 )	/* <Tab> or <Ret> */
	      {	Move_right:
		do { if ( ++icol > nclmn )	icol = 1; }
		while ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 );
		if ( icol == 1 )  goto	Move_down;
		if ( spos != 0 )  goto	Insert;		goto Fillbox;
	      }
	   if ( iret == 12 )				/* case	XK_Up:	  */
	      { Move_up:
		if ( --irow == 0 ) irow = *nrows;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow >= 0)		goto	Move_up;
		if ( spos != 0 )  goto	Insert;	  goto	Fillbox;
	      }
	   if ( iret == 21 )				/* case	XK_Left:  */
	      { Move_left:
		if ( --spos < 0 )
		{ spos = 0;
		     if (--icol == 0 ) { icol = nclmn;   --irow;}
		     if ( irow == 0 )	 icol = irow = 1;
		     if (If_empty(icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 )	goto	Move_left;
		     if ( lvalue != 0 )	goto  Insert;	goto	Fillbox;
		}	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 23 )				/* case	XK_Right: */
	      {	if (spos == 0 && lvalue == 0)	strncpy(stri,array+ind,wbox);
		if (++spos > wbox)	goto	Move_right;
		ixa = ix+spos*wsym;		lvalue = 1;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 32 )				/* case	XK_Down:  */
	      { Move_down:
		if ( ++irow > *nrows)	irow = 1;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow <= *nrows)	goto	Move_down;
		if (spos != 0)	goto Insert;	goto Fillbox;
 	      }
	}
	if ( GetEsc(theEvent.xkey) ) goto EndDialog;
	goto Table_control;
Insert:
	strncpy(array+ind,stri,wbox);   goto Fillbox;

EndDialog:
	RETURNPOINTER_;
	XFlush(theDisplay);
	XDestroyWindow(theDisplay, theWindow);
	XFlush(theDisplay);
	ibcursor = -1;
	mode = 0;
	return 0;
}

/**********************************************************************/
int askgrf_ (title, template, array, len, nrows, ngroup, morow, modex)
			INT_ *len, *nrows, *ngroup, *morow, *modex;
			char title[], template[], array[];
/* The function is called from ASKXGR (file for/surv.f)
            that in turn is called from IFKEY (key "O")

   The same as "asktab", but the 1st column and the column next to "||"
   		are drawn in blue and closed for access
Input:	title	- Title of the table
	template  string defining a structure of the table and its 1st line
		  1st line does not appear if all non-'|' symbols are spaces
	array	- data (numbers or strings) for input and output
	len	- length of "array" element according to description
			in the calling routine (80 in the example below)
	nrows	- number of rows in a table to be created
	ngroup	- if > 0 distance (in rows) between blue separating lines
	morow	- if > 0 separates bottom of the table with a fat blue line
		  this part of the table represents switched off windows,
		  i.e. those with negative number.
   ! Note: this parameter is used for output
	modex	- current radial grid type (input)
*/
{	XEvent	theEvent;
	Window	theWindow;
		int	UpLeftx =2, UpLefty =10,/* window corner location */
		Width =452, Height =480,/* defaut window sizes		*/
		ibox, ibold, hbox,	/* current/old box No. & height	*/
		wsym = 8, hsym = 13,	/* symbol width and height	*/
		wbox, lbox,		/* number of chars in box 	*/
		xshif = 5, yshif = 3,	/* table corner 		*/
		mxfields = 20,	nclmn,	/* max / actual # of columns	*/
		nsta[20],nwid[20],	/* ! <= mxfields=20 are allowed	*/
		nlock[10]={1,0},	/* no more than 10 can be locked*/
		ixa0, iya0, ixa, iya,   /* arrow (textcursor) position	*/
		ind,  inold,	/* address of the selected word in "array" */
		isym, spos=0,	/* abs. and rel. position of symbol in table*/
		lvalue,	ixbox=1,/* 1 if box was changed, 0 otherwise	*/
		iret, ixold, iyold, icold, i, j, nend,
		icol, irow, irow0=0, arrdim, ix, ix1, iy, ihelp ;
	char	stri[132], test='0', str[80];
	char	xlab[6][14] = { {"  \"a\"   [m]  "},{" \"a_N\"  [d/l]"},
				{"\"rho_N\" [d/l]"},{" \"psi\"  [Vs] "},
				{"\"rho_V\" [m]  "},{"\"rho_p\" [d/l]"} };
	strncpy(xlab[0],"  \"a\"   [m]  ",13);
	strncpy(xlab[1]," \"a_N\"  [d/l]",13);
	strncpy(xlab[2],"\"rho_N\" [d/l]",13);
	strncpy(xlab[3]," \"psi\"  [Vs] ",13);
	strncpy(xlab[4],"\"rho_V\" [m]  ",13);
	strncpy(xlab[5],"\"rho_p\" [d/l]",13);
	j = strlen(title);	if ( j > 132 )			{
		printf("%s\n\"%s\"\n%s\n","ASKCOL >>> Title",title,
		     "           is too long");		return (-1);	}
	if  ( *len > 132 )				{
		printf("%s \"%s\" %s\n","ASKCOL >>> Table",title,
		     ". Input string is too long");	return (-1);	}
	nend = strlen(template);	if ( nend > 132 )	{
		printf("%s \"%s\" %s\n","ASKCOL >>> Table",title,
		     ". Requested width is too large");	return (-1);	}
	for ( j=0; j < nend ; j++)
	      if ( template[j] != '|'  &&  template[j] != ' ' ) irow0 = 1;

		/****   Split "template" in blocks  ****/
	for ( j=nclmn=nsta[0]=0; j < nend ; j++ )
	    { if ( template[j] == '|' )
		 { nwid[nclmn] = j-nsta[nclmn];
		   nsta[++nclmn] = j+1;	/*	nlock[nclmn] = 0;
		   if ( template[j+1] == '|' )  nlock[nclmn] = 1;*/
		 }
	      if ( nclmn > mxfields ) {
		printf("%s %d\n","ASKCOL >>>  Too many input fields ",nclmn);
			return (-1);	}
	    }
	nwid[nclmn] = nend-nsta[nclmn];		nclmn++;	arrdim = *len;

	ixa0  = ibold = ihelp = 0;	hbox = hsym+3;	UpLeftx = 2;
	Width = 2*xshif+nend*wsym;	lbox = 3*hbox;
	Height = 2*yshif+(*nrows+irow0)*hbox+lbox;
	if (Kevent.xcur+Kevent.ycur < 10)
		{Kevent.xcur = 330; Kevent.ycur = 430;}
	/*  j = XWX-Width-20;	if (j > UpLeftx)  UpLeftx = j; */
	GetRWgeometry (&XRW,&YRW);
	UpLeftx = 2;	   j = XRW-Width-8;	UpLefty = YRW;
	if (j > UpLeftx)   UpLeftx = j;		j = YRW+Height+30;
	if (j > theHeight) UpLefty = theHeight-Height-30;
	theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0, title, 0,
		    RootWindow(theDisplay,theScreen), theMenuCursor);
	MVPOINTER_ xshif+38*wsym,Height+0.5*hbox-lbox);
	XSelectInput (theDisplay, theWindow, POLL_EV_MASK);
	if ( *modex == 0 )	test = '0';
	if ( *modex == 1 )	test = '1';
	if ( *modex == 2 )	test = '2';
	if ( *modex == 3 )	test = '3';
	if ( *modex == 4 )	test = '4';
	if ( *modex == 5 )	test = '5';

Create_table:
	Change_Color (theGCA, 3, 0); /* white background, blue  foreground */
	if  ( irow0 )
	    {	/*********  Type template in blue  *********/
	    for ( j=0, ix=xshif; j < nclmn ; j++)
		{ WRITE_ theGCA, ix, hsym, template+nsta[j], nwid[j]);
		  if (j+1 == nclmn) break;	ix += wsym*(nwid[j]+1);
		}
		/*********  Draw blue horizontal lines  *********/
	    j = hbox;			LINGCA_ 0,j+4,Width,j+4);
	    }
	if  ( *ngroup > 0 )
	    for ( j = (*ngroup+irow0)*hbox; j+4+yshif < Height-lbox ;
		  j += (*ngroup)*hbox )	  LINGCA_ 0,j+4,Width,j+4);
	if  ( *morow > 0 )
	    {	j = Height-*morow*hbox-lbox-3;	LINGCA_ 0,j,Width,j);
		++j;				LINGCA_ 0,j,Width,j);
	    }
		/*********  Draw blue vertical lines  *********/
	for ( j = 0, ix=xshif-0.5*wsym; j < nclmn-1; j++)
	    { ix += wsym*(nwid[j]+1);	ix1 = ix;
	      if (nwid[j] == 0)		ix1 = ix1-2;
	      if (nwid[j+1] == 0)	ix1 = ix1+2;
	      LINGCA_ ix1,0,ix1,Height-lbox);
	    }
	    j = Height-lbox;	LINGCA_ 0,j,Width,j);
		++j;		LINGCA_ 0,j,Width,j);
	    if ( *modex >= 0 )	LINGCA_ ix1,j,ix1,j+hbox);
		/*********  Draw radial cooordinate box *********/
						/* Change of X-axis */
	   { strcpy(stri,"The current abscissa is ");
	     WRITE_ theGCA, xshif, Height-3+hbox-lbox, stri, 24);
	   }

		/*********  Draw bottom line & comments  *********/
	Change_Color (theGCA, 1, 0);	Change_Color (hghGC, 1, 0);
	WRITE_ theGCA, xshif, Height-3-hbox,
		"Button - select, <Tab>,<Ret>,Arrows - move", 42);
	WRITE_ theGCA, xshif+4*wsym, Height-4, "/<ESC> - done,", 14);
	WRITE_ hghGC,  xshif+28*wsym, Height-4, "?",1);
	WRITE_ theGCA, xshif+29*wsym, Height-4, " - help", 7);
	Change_Color (theGCA, 50, 0); 		/* red on white */
	j = Height-2*hbox; 	LINGCA_ 0,j,Width,j);
	j++; 		   	LINGCA_ 0,j,Width,j);
	Change_Color (hghGC, 50, 0);
	WRITE_ hghGC,  xshif, Height-4, " OK ",4);
	Change_Color (theGCA, 1, 0);	Change_Color (hghGC, 1, 0);
						/* black on white */

		/*********  Draw input data  *********/
	for (   i=0; i < *nrows; i++ )				/* row loop */
 	    {	iy =yshif+hsym+hbox*(i+irow0);	Change_Color(theGCA,3,0);
	    for ( j=isym=0, ix=xshif; j < nclmn ; j++)	     /* column loop */
	    	{ if (j > 0) { isym = nsta[j]; ix += wsym*(nwid[j-1]+1);
			if ( nwid[j-1] == 0)	Change_Color(theGCA,3,0);
			       else		Change_Color(theGCA,1,0);  }
	    	  WRITE_ theGCA, ix, iy, array+i*arrdim+isym, nwid[j]);
		}
	    }
	if ( ibold > 0 && ixbox == 0 )		/* after Expose event */
	   {	isym = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=xshif+wsym*isym;	iy=yshif+hsym+hbox*(irow+irow0-1);
		strncpy(stri,array+ind,wbox);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	Change_Color (hghGC, 3, 0);	ix1 = xshif+24*wsym;
	if ( *modex >= 0 )
	   { WRITE_ hghGC, ix1, Height-3+hbox-lbox, xlab[*modex], 13);
	   }
	else
	   { WRITE_ hghGC, ix1, Height-3+hbox-lbox, "\"time\"", 6);
	   }
	Change_Color(hghGC,1,0);
	if ( ixbox > 0 )		/* Select Xbox on entry */
	   { irow = 1;	ibox = ibold = icol = icold = 0;    lvalue = 0;
	     inold = ind = isym = 0;		wbox = 1;
	     ix = ixold = ixa = ixa0 = xshif+41*wsym-3;
	     iy = iyold = iya = iya0 = Height-3+hbox-lbox;
	     goto	FillXbox;
	   }
	 else if ( *modex >= 0 )
	   { ix1 = xshif+41*wsym;	Change_Color(theGCA,1,0);
	     WRITE_ theGCA, ix1, Height-3+hbox-lbox, &test, 1);
	   }

Table_control:
	XNextEvent (theDisplay, &theEvent);
	if ( theEvent.xany.window == theRootWindow )
{
	     ProcessRootWindowEvent (&theEvent);
	     goto Table_control;
	   }
	if (theEvent.type == Expose)		goto	Create_table;
	if (theEvent.type == ButtonPress)
	{						/* ButtonPressed */
	   ix = theEvent.xbutton.x;	iy = theEvent.xbutton.y;
	   irow = iy-yshif-hbox*irow0-2;	
	   irow = irow/hbox+1;		i=(ix-xshif)/wsym;

	   if ( irow < 0)	goto	Unselect;
	   if ( irow > *nrows+1 && i < 3 && iy > Height-4-hsym )
	      {	if ( spos != 0 )	goto	Escape;	/* OK was pressed */
					goto	EndDialog;
	      }
	   if ( irow == *nrows+1 && i >= 38 && i <= 44 && *modex >= 0 )
	      {	if ( ixbox > 0 )  goto	Table_control;	/* Staying in Xbox */
		ixbox = 1; if ( ibox == 0 ) goto FillXbox;
		goto	Unselect;			/* Entering Xbox */
	      }
NewSelection:						/* Obox selected */
	   for ( j=0; j < nclmn; )
	      { i -= nwid[j]+1;	icol=++j;   if ( i < 0 ) break;	}
	   if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
		   || icol == 1	|| nwid[icol-2] == 0 || nwid[icol-1] == 0 )
	      { ibox = 0;	goto	Unselect;   /* Empty Obox selected */
	      }
	   else
	      { i = (irow-1)*nclmn+icol;	/* Non-empty Obox selected */
		if ( i == ibox)	goto	Table_control;
		if ( i != ibox) { ibold = ibox;	ibox = i;  goto	Unselect; }
	      }
Unselect:
	   if ( ibold > 0 )
		{ strncpy(array+ind,stri,wbox);		ibold = spos = 0;
		  WRITE_ theGCA, ixold, iyold, array+inold, nwid[icold-1]);
		  MoveArrow(theWindow,ixa0,iya0,0,iy);	}
	   else if ( ixbox > 0 )
		{ ixbox = 0;				/* Leaving Xbox */
		  if ( *modex < 0 )	{ ixa0 = 0;    goto	Highlight; }
		  ix1 = xshif+41*wsym;	j = Height-3+hbox-lbox;
		  MoveArrow(theWindow, ixa0, iya0, 0, j);	ixa0 = 0;
		  Change_Color(theGCA,1,0);	WRITE_ theGCA,ix1,j,&test,1);
		  Change_Color (theGCA, 50, 0);	j = Height-2*hbox;
		  LINGCA_ 0,j,Width,j);	j++; 	LINGCA_ 0,j,Width,j);
		  Change_Color(theGCA,1,0);
		}
Highlight:
	   if ( ixbox > 0 )	goto	FillXbox;
	   if ( ibox == 0 )	goto	Table_control;

FillObox:  isym = spos = 0;	if ( icol > 1)  isym = nsta[icol-1];
	   ix=ixa=xshif+wsym*isym;	iy=iya=yshif+hsym+hbox*(irow+irow0-1);
	   ind = (irow-1)*arrdim+isym;	wbox = nwid[icol-1];
	   if ( ibold != 0 )
	        WRITE_ theGCA, ixold, iyold, stri, nwid[icold-1]);
	   if ( ixbox == 0 )
	      { strncpy(stri,array+ind,wbox);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, ixa0, iya0, ixa, iya);
	      }
	   icold = icol;	inold = ind;	ixold = ix;	iyold = iy;
	   ixa0 = ixa;	iya0 = iya;	ibold = ibox;	lvalue = 0;
	   XFlush(theDisplay);
	   goto Table_control;

FillXbox:  if ( *modex < 0 )  goto	Unselect;
	   ix1 = xshif+41*wsym;		iy = Height-3+hbox-lbox;
	   MoveArrow(theWindow, ixa0, iya0, ix1-2, iy);
	   ixa0 = ix1-2;		iya0 = iy;	ibox = 0;
	   sscanf(&test,"%d",modex);	/* Equivalent to atoi(&test) */
	   Change_Color (hghGC, 3, 0);	ix1 = xshif+24*wsym;
	   WRITE_ hghGC, ix1, Height-3+hbox-lbox, xlab[*modex], 13);
 	   ix1 = xshif+41*wsym;	Change_Color(hghGC,1,0);
	   WRITE_ hghGC, ix1, iy, &test, 1);
	   XFlush(theDisplay);
	   goto Table_control;
	}

   if (theEvent.type != KeyPress)	goto	Table_control;
	iret = GetKey (theEvent.xkey, stri, &spos);
	if ( ixbox > 0 )
	   {  					/*  Key pressed in Xbox */
	   if ( iret == 3  || iret == 12 || iret == 23 )
	      { test=test+1;	   if ( test > '5' )	test = '0';
		goto	FillXbox;
	      }
	   if ( iret == 21 || iret == 32 )
	      { test=test-1;	   if ( test < '0' )	test = '5';
		goto	FillXbox;
	      }
	   if ( iret == 0 )
	      {
	      if ( stri[0] == '?' )	goto	Printinfo;
	      if ( stri[0] < '0' )	{
		   stri[0] = '0';	goto	Printinfo;
					}
	      if ( stri[0] > '5')	{
		   stri[0] = '5';	goto	Printinfo;
					}
              strcpy(&test,&stri[0]);
	      goto	FillXbox;
Printinfo:    printf(" The following flux labels are allowed:\n");
	      for (j=0; j<6; j++) printf("       %d   %.13s\n",j,xlab[j]);
	      }
	   }
	if ( ibold > 0 )
	   {  				/*  Key pressed in active Obox */
	   if ( iret == -1 )				/* <Esc> */
Escape:	      { strncpy(array+ind,stri,wbox);   goto EndDialog;
	      }
	   if ( iret == 7 )			   /* <Alt><Esc> */
	      { strncpy(array+ind,stri,wbox);
		*morow=-1;			goto EndDialog;
	      }
	   if ( iret == 0 )
	      {	if (spos == 0 ) { 		/* 1-st entry in the box */
		for (j=1; j < 132; j++) stri[j]=' '; stri[wbox] = '\0';  }
	        if ( stri[spos] == '?' )
		 { printf("%d  %.13s\n",spos,stri);
		printf("\"%.8s\"\n",array+ind);
		printf("\"%s\"\n",str);
		printf("\"%s\"\n",array+ind);
		for(j=0; j<=wbox; j++)
		if (array[ind+j] != ' ') strncat(str,array+ind+j,1);
		printf("\"%s\"\n",str);

		   strcpy (str,"grep -w \"");
		for(j=0; j<=wbox; j++)
		if ( array[ind+j] != ' ' ) strncat(str,array+ind+j,1);
		   strcat (str,"\" tmp/model.txt");
		   printf("%.50s \n",str);
		   if (!ihelp) printf(" No. Scale  Name  Output expression\n");
		   system(str);		ihelp = 1;
		   goto Table_control;
              }
		if (++spos > wbox) spos = wbox;		lvalue = 1;
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	   }
	   if ( iret == 1 )				/* <Del> or <BS> */
	      {	if (--spos >= 0) { strncpy(stri+spos,stri+spos+1,wbox-spos);
				  stri[wbox-1] = ' ';	}
		if (spos < 0) { spos = 0;  strncpy(stri,array+ind,wbox); }
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( (iret == 2 || iret == 3) && ibold > 0 )	/* <Tab> or <Ret> */
	      {	Move_right:
		do { if ( ++icol > nclmn )	icol = 2; }
		while ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
		 	|| nwid[icol-2] == 0 || nwid[icol-1] == 0 );
		if ( icol == 2 )  goto	Move_down;
		if ( spos != 0 )  goto	Insert;   goto FillObox;
	      }
	   if ( iret == 12 )				/* case	XK_Up:	  */
	      { Move_up:
		if ( --irow == 0 ) irow = *nrows;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow >= 0)	goto	Move_up;
		if ( spos != 0 )  goto	Insert;   goto FillObox;
	      }
	   if ( iret == 21 )				/* case	XK_Left:  */
	      { Move_left:
		if ( --spos < 0 )
		{ spos = 0;
		     if (--icol == 1 )	{ icol = nclmn;   --irow; }
		     if ( irow == 0 )	{ irow = 1;	icol = 2; }
		     if ( If_empty(icol, irow, nclmn, nsta, nwid, array, arrdim)
			|| nwid[icol-2]==0 || nwid[icol-1]==0) goto Move_left;
		     if ( lvalue != 0 )	goto  Insert;   goto FillObox;
		}	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 23 )				/* case	XK_Right: */
	      {	if (spos == 0 && lvalue == 0)	strncpy(stri,array+ind,wbox);
		if (++spos > wbox)	goto	Move_right;
		ixa = ix+spos*wsym;		lvalue = 1;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 32 )				/* case	XK_Down:  */
	      { Move_down:
		if ( ++irow > *nrows)	irow = 1;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow <= *nrows)	goto	Move_down;
		if (spos != 0)	goto Insert;	goto	FillObox;
 	      }
	   }
	if ( GetEsc(theEvent.xkey) )		goto	EndDialog;
	goto Table_control;
Insert:
	strncpy(array+ind,stri,wbox);
	goto	FillObox;
EndDialog:
	RETURNPOINTER_;
	XFlush(theDisplay);
	XDestroyWindow(theDisplay, theWindow);
	XFlush(theDisplay);
	ibcursor = -1;
	return (0);
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
        if (k >= 0 && k <= nmamb) j = (int)MAMK[k];
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
int askcol_ (title, template, array, len, nrows, ngroup, morow)
				INT_ *len, *nrows, *ngroup, *morow;
				char title[], template[], array[];
/* The same as "asktab", but the 1st column
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

*/
{	Window	theWindow;
	XEvent	theEvent;
	int	UpLeftx =2, UpLefty =10,/* window corner location	*/
		Width =452, Height =480,/* defaut window sizes		*/
		ibox, ibold, hbox,	/* current/old box No. & height	*/
		wsym = 8, hsym = 13,	/* symbol width and height	*/
		wbox, 			/* number of chars in box 	*/
		xshif = 5, yshif = 3,	/* table corner 		*/
		mxfields = 20,	nclmn,	/* max / actual # of columns	*/
		nsta[20],nwid[20],nend,	/* ! <= mxfields=20 are allowed	*/
		ixa0, iya0, ixa, iya,   /* arrow (textcursor) position	*/
		ind,  inold,	/* address of the selected word in "array" */
		isym, spos=0,	/* abs. and rel. position of symbol in table*/
		lvalue,		/* 1 if box was changed, 0 otherwise	*/
		iret, ixold, iyold, icold, i, j,
		icol, irow, irow0=0, arrdim, ix, ix1, iy, ii=0;
	char	stri[132];
	j = strlen(title);	if ( j > 132 )			{
		printf("%s\n\"%s\"\n%s\n","ASKCOL >>> Title",title,
		     "           is too long");		return 0;	}
	if  ( *len > 132 )				{
		printf("%s \"%s\" %s\n","ASKCOL >>> Table",title,
		     ". Input string is too long");	return 0;	}
	nend = strlen(template);	if ( nend > 132 )	{
		printf("%s \"%s\" %s\n","ASKCOL >>> Table",title,
		     ". Requested width is too large");	return 0;	}
	for ( j=0; j < nend ; j++)
	      if ( template[j] != '|'  &&  template[j] != ' ' ) irow0 = 1;

		/****   Split "template" in blocks  ****/
	for ( j=nclmn=nsta[0]=0; j < nend ; j++ )
	    { if ( template[j] == '|' )
		 { nwid[nclmn] = j-nsta[nclmn]; nsta[++nclmn] = j+1;}
	      if ( nclmn > mxfields ) {
		printf("%s %d\n","ASKCOL >>>  Too many input fields ",nclmn);
			return 0;	}
	    }
	nwid[nclmn] = nend-nsta[nclmn];		nclmn++;	arrdim = *len;

	ixa0  = ibold = 0;		hbox   = hsym+3;	UpLeftx = 2;
	Width = 2*xshif+nend*wsym;	Height = 2*yshif+(*nrows+2+irow0)*hbox;
	if (Kevent.xcur+Kevent.ycur < 10)
			{Kevent.xcur = 330; Kevent.ycur = 430;}
	/*  j = XWX-Width-20;	if (j > UpLeftx)  UpLeftx = j; */
	GetRWgeometry (&XRW,&YRW);
	UpLeftx = 2;	   j = XRW-Width-8;	UpLefty = YRW;
	if (j > UpLeftx)   UpLeftx = j;		j = YRW+Height+30;
	if (j > theHeight) UpLefty = theHeight-Height-30;
	theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0, title, 0,
		    RootWindow(theDisplay,theScreen), theMenuCursor);
	MVPOINTER_ xshif+wsym*nwid[0]/2,yshif+irow0*hbox+hsym);
	XSelectInput (theDisplay, theWindow, POLL_EV_MASK);

Create_table:
	Change_Color (theGCA, 3, 0); /* white background, blue  foreground */
	if  ( irow0 )
	    {	/*********  Type template in blue  *********/
	    for ( j=0, ix=xshif; j < nclmn ; j++)
		{ WRITE_ theGCA, ix, hsym, template+nsta[j], nwid[j]);
		  if (j+1 == nclmn) break;	ix += wsym*(nwid[j]+1);
		}
		/*********  Draw blue horizontal lines  *********/
	    j = hbox;			LINGCA_ 0,j+4,Width,j+4);
	    }
	if  ( *ngroup > 0 )
	    for ( j = (*ngroup+irow0)*hbox; j+4+yshif < Height-2*hbox ;
		  j += (*ngroup)*hbox )	  LINGCA_ 0,j+4,Width,j+4);
	if  ( *morow > 0 )
	    {	j = Height-(*morow+2)*hbox-3;	LINGCA_ 0,j,Width,j);
		++j;				LINGCA_ 0,j,Width,j);
	    }
		/*********  Draw blue vertical lines  *********/
	for ( j = 0, ix=xshif-0.5*wsym; j < nclmn-1; j++)
	    { ix += wsym*(nwid[j]+1);	ix1 = ix;
	      if (nwid[j] == 0)		ix1 = ix1-2;
	      if (nwid[j+1] == 0)	ix1 = ix1+2;
	      LINGCA_ ix1,0,ix1,Height-2*hbox);
	    }
		/*********  Draw bottom line & comments  *********/
	WRITE_ theGCA, xshif, Height-3-hbox,
		"Button - select, <Tab>,<Ret>,Arrows - move", 42);
	WRITE_ theGCA, xshif+4*wsym, Height-4, "/<ESC> - done", 13);
	Change_Color (theGCA, 50, 0);  /* white background, red foreground */
	j = Height-2*hbox; 	LINGCA_ 0,j,Width,j);
	j++; 		   	LINGCA_ 0,j,Width,j);
	Change_Color (hghGC, 50, 0);
	WRITE_ hghGC,  xshif, Height-4, " OK ",4);
	Change_Color (theGCA, 1, 0);	Change_Color (hghGC, 1, 0);
				/* white background, black foreground */
		/*********  Draw input data  *********/
	for (   i=0; i < *nrows; i++ )				/* row loop */
 	    {	iy =yshif+hsym+hbox*(i+irow0);	Change_Color (theGCA, 3, 0);
	    for ( j=isym=0, ix=xshif; j < nclmn ; j++)	     /* column loop */
	    	{ if (j > 0) { isym = nsta[j]; ix += wsym*(nwid[j-1]+1); }
	    	  WRITE_ theGCA, ix, iy, array+i*arrdim+isym, nwid[j]);
		  Change_Color (theGCA, 1, 0);
		}
	    }
	if ( ibold > 0 )			/* after Expose event */
	   {	isym = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=xshif+wsym*isym;	iy=yshif+hsym+hbox*(irow+irow0-1);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	if ( ii == 0 )		/* Select Upper Left box on entry */
	   { irow =ii = 1;	ibox = ibold = icol = icold = 2;    lvalue = 0;
		inold = ind = isym = nsta[1];		wbox = nwid[1];
		ix = ixold = ixa = ixa0 = xshif+wsym*isym;
		iy = iyold = iya = iya0 = yshif+hsym+hbox*irow0;
		WRITE_ hghGC, ix, iy, array+ind, wbox);
		strncpy(stri,array+ind,wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	XFlush(theDisplay);

Table_control:
	XNextEvent (theDisplay, &theEvent);
	if ( theEvent.xany.window == theRootWindow )
           { ProcessRootWindowEvent (&theEvent);
	     goto Table_control;
	   }
	if (theEvent.type == Expose)	goto	Create_table;
	if (theEvent.type == ButtonPress)
	{   ix = theEvent.xbutton.x;	iy = theEvent.xbutton.y;
	    irow = iy-yshif-hbox*irow0-2;
	    if ( irow >= 0)
	    {	irow = irow/hbox+1;		i=(ix-xshif)/wsym;
		if ( ibold > 0 )	strncpy(array+ind,stri,wbox);
		if ( irow > *nrows )		  /* No selection or Exit */
		   { if ( i > 3 || iy < Height-4-hsym ) goto Unselect;
			goto EndDialog;			/* OK was pressed */
		   }
		for ( j=0; j < nclmn; )			/* Selection is made */
		    { i -= nwid[j]+1;	icol=++j;   if ( i < 0 ) break;	}
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
		   || icol == 1	|| nwid[icol-1] == 0 )	goto	Unselect;
Fillbox:	isym = spos = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=ixa=xshif+wsym*isym;	iy=iya=yshif+hsym+hbox*(irow+irow0-1);
		ind = (irow-1)*arrdim+isym;	wbox = nwid[icol-1];
		ibox = (irow-1)*nclmn+icol;
		if ( ibold != 0 )
		     WRITE_ theGCA, ixold, iyold, stri, nwid[icold-1]);
		strncpy(stri,array+ind,wbox);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, ixa0, iya0, ixa, iya);
		icold = icol;	inold = ind;	ixold = ix;	iyold = iy;
		ixa0 = ixa;	iya0 = iya;	ibold = ibox;	lvalue = 0;
	    }
	    else
Unselect:   {	if ( ibold > 0 )			{
		WRITE_ theGCA, ixold, iyold, array+inold, nwid[icold-1]);
		MoveArrow(theWindow,ixa0,iya0,0,iy);	ibold=spos=0;	}
	    }
	    XFlush(theDisplay);
	    goto Table_control;
	}
	if ( ibold > 0 )
	{
	   iret = GetKey (theEvent.xkey, stri, &spos);
	   if ( iret == -1 && ibold > 0 )			/* <Esc> */
	      { strncpy(array+ind,stri,wbox);   goto EndDialog;
	      }
	   if ( iret == 0 )
	      {	if (spos == 0 ) { 		/* 1-st entry in the box */
		for (j=1; j < 132; j++) stri[j]=' '; stri[wbox] = '\0';  }
		if (++spos > wbox) spos = wbox;		lvalue = 1;
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 1 )				/* <Del> or <BS> */
	      {	if (--spos >= 0) { strncpy(stri+spos,stri+spos+1,wbox-spos);
				  stri[wbox-1] = ' ';	}
		if (spos < 0) { spos = 0;  strncpy(stri,array+ind,wbox); }
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( (iret == 2 || iret == 3) && ibold > 0 )	/* <Tab> or <Ret> */
	      {	Move_right:
		do { if ( ++icol > nclmn )	icol = 2; }
		while ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 );
		if ( icol == 2 )  goto	Move_down;
		if ( spos != 0 )  goto	Insert;		goto Fillbox;
	      }
	   if ( iret == 12 )				/* case	XK_Up:	  */
	      { Move_up:
		if ( --irow == 0 ) irow = *nrows;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow >= 0)		goto	Move_up;
		if ( spos != 0 )  goto	Insert;	  goto	Fillbox;
	      }
	   if ( iret == 21 )				/* case	XK_Left:  */
	      { Move_left:
		if ( --spos < 0 )
		{ spos = 0;
		     if (--icol == 1 )	{ icol = nclmn;   --irow; }
		     if ( irow == 0 )	{ irow = 1;	icol = 2; }
		     if (If_empty(icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 )	goto	Move_left;
		     if ( lvalue != 0 )	goto  Insert;	goto	Fillbox;
		}	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 23 )				/* case	XK_Right: */
	      {	if (spos == 0 && lvalue == 0)	strncpy(stri,array+ind,wbox);
		if (++spos > wbox)	goto	Move_right;
		ixa = ix+spos*wsym;		lvalue = 1;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 32 )				/* case	XK_Down:  */
	      { Move_down:
		if ( ++irow > *nrows)	irow = 1;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow <= *nrows)	goto	Move_down;
		if (spos != 0)	goto Insert;	goto Fillbox;
 	      }
	}
	if ( GetEsc(theEvent.xkey) ) goto EndDialog;
	goto Table_control;
Insert:
	strncpy(array+ind,stri,wbox);   goto Fillbox;

EndDialog:
	RETURNPOINTER_;
	XFlush(theDisplay);
	XDestroyWindow(theDisplay, theWindow);
	XFlush(theDisplay);
	ibcursor = -1;
	return (ibold);
}

/**********************************************************************/
int If_empty(int icol, int irow, int nclmn, int nsta[], int nwid[], char array[], int arrdim){
/* True if cursor is in empty box
  icol - current column
  irow - current row
  nclmn - total # of columns
*/

    int j, imin, imax;
    imin = 1;
    imax = nclmn;

    R1: for (j=imin-1; j<nclmn ; j++){
            if (nwid[j] == 0){
                if (j >= icol){
                    imax = j;
                    break;
                }
                imin = j + 2;
                goto R1;
            }
        }

    if (strspn(array+(irow-1)*arrdim+nsta[imin-1], " ") > nsta[imax-1]+nwid[imax-1]-nsta[imin-1]){
        return 1;
    }
    else{
        return 0;
    }
}

/**********************************************************************/
int GetValue(XKeyEvent theEvent, char str[], int lvalue, char tstri[], int *pos){
    int i;
    XComposeStatus theComposeStatus;
    KeySym  theKeySym;
    int theBufferLength, theKeyBufferMaxLen = 4;
    char theKeyBuffer[5];

    theBufferLength = XLookupString(&theEvent, theKeyBuffer,
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
     case XK_Return:
     case XK_KP_Enter: return 1;
     case XK_Tab:      return 2;
     case XK_Escape:   return -1;
     case XK_question: return -2;
     case XK_R7:       return 11;
     case XK_Left:     return 21;
     case XK_Up:       return 12;
     case XK_Right:    return 23;
     case XK_Down:     return 32;
     case XK_R9:       return 13;
     case XK_R13:      return 31;
     case XK_R15:      return 33;
     default:
         if (*pos<lvalue){
             for (i=0; tstri[i] != 0; i++){
                 if (theKeySym == tstri[i]){
                     if (*pos < 0) *pos=0;
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
int GetKey(XKeyEvent theKeyEvent, char str[], int *pos){
    XComposeStatus theComposeStatus;
    KeySym theKeySym;
    int theKeyBufferMaxLen = 4;
    char theKeyBuffer[5];

    XLookupString(&theKeyEvent, theKeyBuffer,
        theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
    if (theKeyEvent.state & Mod1Mask){
        if (theKeySym == XK_Escape){
            return 7;
        }
    }

    switch(theKeySym){
    case XK_Escape:    return -1;
    case XK_BackSpace:
    case XK_Delete:    return 1;
    case XK_Return:
    case XK_KP_Enter:  return 2;
    case XK_Tab:       return 3;
    case XK_R7:        return 11;
    case XK_Left:      return 21;
    case XK_Up:        return 12;
    case XK_Right:     return 23;
    case XK_Down:      return 32;
    case XK_R9:        return 13;
    case XK_R13:       return 31;
    case XK_R15:       return 33;
    case XK_Shift_L:               /* Left shift */
    case XK_Shift_R:   return 34;  /* Right shift */
    case XK_Control_L:             /* Left control */
    case XK_Control_R: return 35;  /* Right control*/
    case XK_Meta_L:    break;      /* Left meta */
    case XK_Meta_R:    return 36;  /* Right meta */
    case XK_Alt_L:                 /* Left alt */
    case XK_Alt_R:     return 37;  /* Right alt */
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
int GetEsc(XKeyEvent theKeyEvent){
    XComposeStatus theComposeStatus;
    KeySym theKeySym;
    int theKeyBufferMaxLen=4;
    char theKeyBuffer[5];

    XLookupString(&theKeyEvent, theKeyBuffer, theKeyBufferMaxLen,
        &theKeySym, &theComposeStatus);
    if (theKeySym == XK_Escape) return 1;
    return 0;
}

/**********************************************************************/
int GetName(XKeyEvent theKeyEvent, char str[], int lvalue, int *pos){
    XComposeStatus theComposeStatus;
    KeySym theKeySym;
    int theKeyBufferMaxLen = 4, i;
    char theKeyBuffer[5], lsym;

    XLookupString(&theKeyEvent, theKeyBuffer,
        theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
    switch(theKeySym){
    case XK_BackSpace:
    case XK_Delete:
        *pos -= 1;
        if (*pos >= 0) str[*pos] = ' ';
        break;
    case XK_Return:
    case XK_KP_Enter: return 1;
    case XK_Tab:      return 2;
    case XK_Escape:   return -1;
    case XK_R7:       return 11;
    case XK_Left:     return 21;
    case XK_Up:       return 12;
    case XK_Right:    return 23;
    case XK_Down:     return 32;
    case XK_R9:       return 13;
    case XK_R13:      return 31;
    case XK_R15:      return 33;
    default:
        if (theKeySym > 127) return 0;
        if (*pos < lvalue){
            if (*pos < 0) *pos = 0;     /* Restore default name */
            for (i=*pos; i<lvalue; i++) str[i]=' ';
            lsym = theKeySym;
            if (isalnum(lsym) || lsym == '/' || lsym == '.' || lsym == '_'){
                str[*pos] = lsym;
                *pos += 1;
            }
        }
        break;
    }
    return 0;
}

/**********************************************************************/
int ufilebox_(INT_ *nofbox, char una[], char theNames[], char unad[]){
/* U-file name setting */
    Window theWindow;
    XEvent theEvent;
    int UpLeftx=2, UpLefty=10,  /* window corner location       */
        Width=452, Height=480;  /* window size                  */
    int ix, iy, wbox, hbox,     /* parameter box corner and size*/
        wsym=8, hsym=13,        /* symbol width and height      */
        lvalue=10, lname=4,     /* length of value and name     */
        xshif=5, yshif=3,       /* table corner and             */
        nclmn=3, nparam,        /* # of columns and parameters  */
        ixa0=0, iya0=0, ixa,iya,/* arrow (textcursor) position  */
        i, j, ibox, spos = 0, vpos, lsym=4,
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
    nparam = *nofbox;
    Width = 2*xshif + nclmn*wbox - wsym;
    Height =2*yshif + ((nparam - 1)/nclmn + 3)*hbox;
    ibox = 0;
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
        ix = wbox*((i - 1)%nclmn) + xshif;
        iy = hbox*((i - 1)/nclmn) + yshif;
        ind = (i - 1)*lname;
        WRITE_ theGCA, ix, hsym+iy, theNames+ind, lname);
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
    for (i=1; i < nclmn; i++){
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
    i = nclmn*wbox/wsym-5;
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
            if (ibox){
                oldparam = ibox;
                ixold = ix;
                iyold = iy;
            }
            ibox = FindBoxNum(xButton-xshif+wsym/2, yButton-yshif, wbox, hbox, nclmn, nparam);
            i = FindBoxNum(xButton-xshif, yButton-Height+hbox+1, 4*wsym, hbox, 1, 1);
            if (i){
                iret = -1;
                ibox = oldparam;
                ixold = ix;
                iyold = iy;
                key_flag = 1;
            }
            else{
                if (ibox == 0) continue;
                ind = (oldparam - 1)*lvalue;
                i = ibox;
                ix = wbox*((i - 1)%nclmn) + xshif;
                iy = hbox*((i - 1)/nclmn) + yshif;
                ind = (oldparam - 1)*lname;
                ind1 = (i - 1)*lname;
                WRITE_ theGCA, ixold, hsym+iyold, theNames+ind, lname);
                WRITE_ hghGC, ix, hsym+iy, theNames+ind1, lname);
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
            if (ibox == 0){
                if (GetEsc(theEvent.xkey)){
                    break; // destroy u-file menu window
                }
                else{
                    continue;
                }
            }
            else{
                if (theEvent.type == KeyPress){
                    ind = (ibox-1)*lvalue;
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
            i = nclmn*wbox/wsym - 5;
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
            oldparam = ibox;
            if (iret ==  2){
                ibox++;
                if (ibox > nparam) ibox = 1;
            }
            else{
                stcopy(Ufile_Name+4, una+ind, lvalue);
                dcol = iret%10 - 2;
                icol = (ibox - 1)%nclmn + dcol;
                if (icol < 0) icol = 0;
                if (icol >= nclmn) icol--;
                drow = iret/10 - 2;
                irow = (ibox - 1)/nclmn + drow;
                if (irow < 0) irow = 0;
                if (irow > (nparam-1)/nclmn) irow--;
                ibox = irow*nclmn + icol + 1;
                if (ibox == oldparam && icol == nclmn-1 && drow != -1) ibox++;
                if (ibox == oldparam && icol == 0       && drow !=  1) ibox--;
                if (ibox < 1) ibox = 1;
                if (ibox > nparam) ibox = nparam;
            }
            ixold = ix;
            iyold = iy;
            spos = 0;
            i = ibox;
            ix = wbox*((i - 1)%nclmn) + xshif;
            iy = hbox*((i - 1)/nclmn) + yshif;
            ind = (oldparam - 1)*lname;
            ind1 = (i - 1)*lname;
            WRITE_ theGCA, ixold, hsym+iyold, theNames+ind, lname);
            WRITE_ hghGC, ix, hsym+iy, theNames+ind1, lname);
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

    RETURNPOINTER_;
    XFlush(theDisplay);
    XDestroyWindow(theDisplay, theWindow);
    ibcursor = -1;
    XFlush(theDisplay);
    return 0;
}

/**********************************************************************/
int FindBoxNum(int ix, int iy, int wBox, int hBox, int nclmn, int nBox){
/* Returns box # or 0
  ix, iy:          current coordinates
  nclmn, nBox:     #columns, #boxes
  wBox, hBox:      box width, height
*/
    int icol, pnum;

    if (ix <= 0 || iy < 0) return 0;
    icol = ix/wBox;
    if (icol >= nclmn) return 0;
    pnum = (iy/hBox)*nclmn + icol + 1;
    if (pnum > nBox) return 0;
    return pnum;
}
