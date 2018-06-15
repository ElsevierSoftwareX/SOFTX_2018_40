%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  KINEMATIC COLLAPSE LOAD CALCULATOR                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB script in combination with the MATLAB figure
% 'KCLCalculator.fig' produce a GUI for the static analysis of
% circular masonry arches subjected to a point load or constant horizontal
% acceleration condition. 
% The GUI was developed by:
%                           Gabriel Stockdale
%                           Simone Tiberti, 
%                           Daniela Camilletti,
%                           Gessica Sferrazza Papa
%                           Ahmad Basshofi Habieb
%                           Elisa Bertolesi
%                           Gabriele Milani
%                           Siro Casolo
%                                                           March 13, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY
% The following script and accompanying figure file result in a GUI that 
% creates an interactive platform for the examination of the collapse load 
% and hinge configuration of a circular masonry arch for an asymmetric 
% point load and horizontal acceleration conditions.
% The collapse load and hinge reactions are calculated through the use of
% an equilibrium method to the kinematic approach of limit analysis.
% After running the script, the user inputs the arch data including the
% number of blocks, internal radius, thickness to radius ratio, depth,
% density, and load case. After inputting these values and pressing the run
% button, the program establishes the arch geometry and weight, sets the
% four hinge locations, calculates the collapse load and hinge reactions,
% and then plots the thrust line. The user then adjusts the hinge sliders
% to change the mechanism. Each time the mechanism changes the calculations
% are performed and the results are plotted. The thrust line is calculated
% by determining the necessary eccentricity from the centerline of the arch
% that is required at each point such that the moment is zero when 
% evaluated at hinge 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- INITIALAZING GUI -----------------------------------
function varargout = masonry_settembre(varargin)
% KCLCALCULATOR MATLAB code for KCLCalculator.fig
%      KCLCALCULATOR, by itself, creates a new 
%           KCLCALCULATOR or raises the existing singleton*.
%
%      H = KCLCALCULATOR returns the handle to a new 
%          KCLCALCULATOR or the handle to the existing singleton*.
%
%      KCLCALCULATOR('CALLBACK',hObject,eventData,handles,...)
%       calls the local function named CALLBACK in 
%       KCLCALCULATOR.M with the given input arguments.
%
%      KCLCALCULATOR('Property','Value',...) creates a new 
%       KCLCALCULATOR or raises the existing singleton*.  
%       Starting from the left, property value pairs are applied to the GUI
%       before KCLCalculator_OpeningFcn gets called.  An
%       unrecognized property name or invalid value makes property 
%       application stop.  All inputs are passed to 
%       KCLCalculator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
 
% Edit the above text to modify the response to help KCLCalculator
 
% Last Modified by GUIDE v2.5 13-Mar-2018 16:46:10
 
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KCLCalculator_OpeningFcn,...
                   'gui_OutputFcn',  @KCLCalculator_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
 
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before KCLCalculator is made visible.
function KCLCalculator_OpeningFcn(hObject, eventdata,...
    handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure_1
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KCLCalculator(see VARARGIN)
 
 
% Choose default command line output for KCLCalculator
handles.output = hObject;
 
% Update handles structure
guidata(hObject, handles);
 
%%aggiungo io
set(hObject,'toolbar','figure'); %per modificare le figure
 
% UIWAIT makes KCLCalculator wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout =...
    KCLCalculator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure_1
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Get default command line output from handles structure
varargout{1} = handles.output;
 
% -------------------------------------------------------------------------
%% ---------------------- HINGE SLIDERS -----------------------------------
% The hinge sliders change the location of their respective hinge. The
% hinge limits are:
% Hinges H1 and H4 are the primary hinges and their limits are the base and
% two joints below the crown on their respective sides.
% Hinges H2 and H3 are each bound between H1 and H4 respectively and the
% crown.
% The Eval function is called at the end of each slider callback function.
%  ------------------------ Slider H1 -------------------------------------
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain data from handles ------------------------------------------------
% N - Number of blocks
N = str2num(get(handles.editText_N, 'String'));
% R - Internal radius
R = str2num(get(handles.editText_R,'String')); % [m]
% T_R - Thickness to radius ratio
T_R = str2num(get(handles.editText_T_R,'String')); 
% Hinge H2 is called for comparison
H2 = get(handles.slider2,'Value');
% H1max - Max value of slider for evaluations
H1max = round(get(handles.slider1,'Max'));
% -------------------------------------------------------------------------
% Base calculations -------------------------------------------------------
% Thickness
T = T_R*R; %[m]
% Extrados Radius
RE = R+T; %[m]
% Bndry_blocks - The step between block boundaries
Bndry_blocks = pi/N;
% -------------------------------------------------------------------------
% Read the H1 slider value ------------------------------------------------
H1 = round(get(hObject, 'Value'));
% -------------------------------------------------------------------------
% Calculate x and y coordinates of H1 -------------------------------------
H1x = RE*cos(Bndry_blocks*H1); %[m]
H1y = RE*sin(Bndry_blocks*H1); %[m]
% -------------------------------------------------------------------------
% Check hinge 1 position against hinge 2 ----------------------------------
% This keeps H2 ahead of H1 always
if H1 >= H2
    H2 = H1+1;
    H2x = R*cos(Bndry_blocks*H2); %[m]
    H2y = R*sin(Bndry_blocks*H2); %[m]
    set(handles.marker2,'Xdata',H2x,'Ydata',H2y)
    set(handles.marker2txt,'Position',[H2x H2y])
end
% H2 slider disabled when H1 is at its maximum
if H1 == H1max
    set(handles.slider2,'Enable','off')
else
    set(handles.slider2,'Enable','on')
end
% Update H2 and its lower limit
set(handles.slider2,'Value',H2)
set(handles.slider2,'Min',H1+1)
% Manage sliderstep for H2. This removed a locking of the slider that was
% observed when H2 had two possible locations. This forces the sliderstep
% to have two choices at that time.
if H1 == H1max-1
    set(handles.slider2,'SliderStep',[1 1])
else
    set(handles.slider2,'SliderStep',[(1/((N-1)/2-H1)) 1])
end
% -------------------------------------------------------------------------
% Update the marker and slider value for H1 -------------------------------
set(handles.marker1,'Xdata',H1x,'Ydata',H1y)
set(handles.marker1txt,'Position',[H1x H1y])
set(hObject, 'Value', H1);
% -------------------------------------------------------------------------
% Run the Eval function ---------------------------------------------------
Eval(hObject, eventdata, handles)
% -------------------------------------------------------------------------
% End slider1_Callback ----------------------------------------------------
  
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% -------------------------------------------------------------------------
% End slider1_CreatFcn ----------------------------------------------------
% ----------------------- End Slider H1 -----------------------------------
 
% ------------------------- Slider H2 ------------------------------------
% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Obtain data from handles ------------------------------------------------
% N - Number of blocks
N = str2num(get(handles.editText_N, 'String'));
% R - Internal radius
R = str2num(get(handles.editText_R,'String')); % [m]
% T_R - Thickness to radius ratio
T_R = str2num(get(handles.editText_T_R,'String')); 
% -------------------------------------------------------------------------
% Base calculations -------------------------------------------------------
% Thickness
T = T_R*R; % [m]
% Extrados Radius
RE = R+T; % [m]
% Bndry_blocks - The step between block boundaries
Bndry_blocks = chop(pi/N,4);
% -------------------------------------------------------------------------
% Read the H2 slider value ------------------------------------------------
H2 = round(get(hObject, 'Value'));
% -------------------------------------------------------------------------
% Calculate x and y coordinates of H2 -------------------------------------
H2x = R*cos(Bndry_blocks*H2); %[m]
H2y = R*sin(Bndry_blocks*H2); %[m]
% -------------------------------------------------------------------------
% Update the marker and slider value for H2 -------------------------------
set(handles.marker2,'Xdata',H2x,'Ydata',H2y)
set(handles.marker2txt,'Position',[H2x H2y])
set(hObject, 'Value', H2);
% -------------------------------------------------------------------------
% Run the Eval function ---------------------------------------------------
Eval(hObject, eventdata, handles)
% -------------------------------------------------------------------------
% End slider2_Callback ----------------------------------------------------
    
% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
if isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end   
% -------------------------------------------------------------------------
% End slider2_CreatFcn ----------------------------------------------------
% ----------------------- End Slider H2 -----------------------------------
 
% ------------------------- Slider H3 -------------------------------------
% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain data from handles ------------------------------------------------
% N - Number of blocks
N = str2num(get(handles.editText_N, 'String'));
% R - Internal radius
R = str2num(get(handles.editText_R,'String')); % [m]
% T_R - Thickness to radius ratio
T_R = str2num(get(handles.editText_T_R,'String'));
% Get the point load checkbox value
Pl=get(handles.checkbox_PL,'Value');
% -------------------------------------------------------------------------
% Base calculations -------------------------------------------------------
% Thickness
T = T_R*R; % [m]
% Extrados Radius
RE = R+T; % [m]
% Bndry_blocks - The step between block boundaries
Bndry_blocks = pi/N;
% -------------------------------------------------------------------------
% Read the H3 slider value ------------------------------------------------
H3 = round(get(hObject, 'Value'));
% -------------------------------------------------------------------------
% Calculate x and y coordinates of H3 -------------------------------------
H3x = RE*cos(Bndry_blocks*H3); %[m]
H3y = RE*sin(Bndry_blocks*H3); %[m]
% -------------------------------------------------------------------------
% Update the marker and slider value for H3 -------------------------------
set(handles.marker3,'Xdata',H3x,'Ydata',H3y)
set(handles.marker3txt,'Position',[H3x H3y])
set(hObject, 'Value', H3);
% Include the point load arrow with marked
if Pl == 1
    set(handles.PLt,'Position',[H3x H3y])
end
% -------------------------------------------------------------------------
% Run the Eval function ---------------------------------------------------
Eval(hObject, eventdata, handles)
% -------------------------------------------------------------------------
% End slider3_Callback ----------------------------------------------------
 
% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if isequal(get(hObject,'BackgroundColor'),...
                                get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% -------------------------------------------------------------------------
% End slider3_CreatFcn ----------------------------------------------------
% ----------------------- End Slider H3 -----------------------------------
 
% ------------------------- Slider H4 -------------------------------------
% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Obtain data from handles ------------------------------------------------
% N - Number of blocks
N = str2num(get(handles.editText_N, 'String'));
% R - Internal radius
R = str2num(get(handles.editText_R,'String')); % [m]
% T_R - Thickness to radius ratio
T_R = str2num(get(handles.editText_T_R,'String')); 
% Hinge H2 is called for comparison
H3 = get(handles.slider3,'Value');
% H1max - Max value of slider for evaluations
H4min = round(get(handles.slider4,'Min'));
% Get the point load checkbox value
Pl=get(handles.checkbox_PL,'Value');
% -------------------------------------------------------------------------
% Base calculations -------------------------------------------------------
% Thickness
T = T_R*R; %[m]
% Extrados Radius
RE = R+T; %[m]
% Bndry_blocks - The step between block boundaries
Bndry_blocks = pi/N;
% -------------------------------------------------------------------------
% Read the H4 slider value ------------------------------------------------
H4=round(get(hObject, 'Value'));
% -------------------------------------------------------------------------
% Calculate x and y coordinates of H4 -------------------------------------
H4x = R*cos(Bndry_blocks*H4); %[m]
H4y = R*sin(Bndry_blocks*H4); %[m]
% -------------------------------------------------------------------------
% Check H4 position against H3 --------------------------------------------
% This keeps H4 ahead of H3 always (polar)
if H4 <= H3
    H3 = H4-1;
    H3x = RE*cos(Bndry_blocks*H3); %[m]
    H3y = RE*sin(Bndry_blocks*H3); %[m]
    set(handles.marker3,'Xdata',H3x,'Ydata',H3y)
    set(handles.marker3txt,'Position',[H3x H3y])
    if Pl ==1
        set(handles.PLt,'Position',[H3x H3y])
    end
end
% H3 slider disabled when H4 is at its minimum
if H4 == H4min
    set(handles.slider3,'Enable','off')
else
    set(handles.slider3,'Enable','on')
end
% Update H3 and its upper limit
set(handles.slider3,'Max',H4-1)
set(handles.slider3,'Value',H3)
% Manage sliderstep for H3. This removed a locking of the slider that was
% observed when H2 had two possible locations. This forces the sliderstep
% to have two choices at that time.
if H4 == H4min+1
    set(handles.slider3,'SliderStep',[1 1])
else
    set(handles.slider3,'SliderStep',[(1/(H4-((N+1)/2))) 1])
end 
% -------------------------------------------------------------------------
% Update the marker and slider value for H1 -------------------------------
set(handles.marker4,'Xdata',H4x,'Ydata',H4y)
set(handles.marker4txt,'Position',[H4x H4y])
set(hObject, 'Value', H4);
% -------------------------------------------------------------------------
% Run the Eval function ---------------------------------------------------
Eval(hObject, eventdata, handles)
% -------------------------------------------------------------------------
% End slider4_Callback ----------------------------------------------------
 
% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% -------------------------------------------------------------------------
% End slider4_CreatFcn ----------------------------------------------------
% ----------------------- End Slider H4 -----------------------------------
% --------------------- END HINGE SLIDERS ---------------------------------
 
 
%% ----------------------- RUN BUTTON -------------------------------------
% The run push button performs three functions. First, the input parameters
% are checked to make sure they exist and fall within the prescribed
% limits. Second, the arch is drawn, the initial hinge locations are set
% and marked on the figure. Third, the Eval function is called to perform
% the analysis calculations, display the load and reaction values, and plot
% the thrust line.
% --- Executes on button press in run_pushbutton.
function run_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to run_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Clear the figure --------------------------------------------------------
cla(handles.figure_1, 'reset')
% -------------------------------------------------------------------------
% Obtain data from handles ------------------------------------------------
% N - Number of blocks
N = str2num(get(handles.editText_N, 'String'));
% R - Internal radius
R = str2num(get(handles.editText_R,'String')); % [m]
% T_R - Thickness to radius ratio
T_R = str2num(get(handles.editText_T_R,'String')); 
% Block depth
Depth_blocks=str2num(get(handles.editText_Depth,'String')); %[m]
% Block density
Den_blocks=str2num(get(handles.editText_Density,'String')); %[kg/m^3]
% Hinge H3 is called for comparison
H3 = get(handles.slider3,'Value');
% H4min - Min value of slider for evaluations
H4min = round(get(handles.slider4,'Min'));
% Point load checkbox value
Pl=get(handles.checkbox_PL,'Value');
% Horizontal acceleration checkbox value
H_acc=get(handles.checkbox_H_acc,'Value');
% -------------------------------------------------------------------------
% Base calculations -------------------------------------------------------
% Thickness
T = T_R*R; %[m]
% Extrados Radius
RE = R+T; %[m]
% Bndry_blocks - The step between block boundaries
Bndry_blocks = pi/N;
% -------------------------------------------------------------------------
% Set load type -----------------------------------------------------------
if Pl == 1
    loadType = 1;
elseif H_acc == 1
    loadType = 2;
end 
% -------------------------------------------------------------------------
% Check for valid input data ----------------------------------------------
% Requirements:
%   i - It must exist
%  ii - The number of blocks must be odd and greater than 5
% iii - Thickness to radius ratio must be between 0.11 and 0.33
check_N= mod(N,2); 
if isempty(N)||isempty(R)||isempty(T_R)||isempty(Depth_blocks)...
        ||isempty(Den_blocks)
   errordlg('Input data is missing','Error');
elseif check_N==0
   errordlg('N must be odd','Error');
elseif N<5
   errordlg('N must be greater than 5','Error');
elseif T_R<=0.11
    errordlg('T_R must be greater than 0.11','Error');
elseif T_R>=0.33
    errordlg('T_R must be less than 0.33','Error');
elseif (Pl==1 && H_acc==1)||(Pl==0 && H_acc==0)
    errordlg('Select a load type', 'Error');
else 
% valid input data --------------------------------------------------------
% -------------------------------------------------------------------------
% ---------------------------- PLOT FIGURE --------------------------------
    % Set figure 1 to current axes and set parameters ---------------------
    axes(handles.figure_1);
    % Fix the aspect ratio
    set(gca,'DataAspectRatio',[1 1 1])
    % Set axis labels
    xlabel('Length [m]');
    ylabel('Height [m]');
    hold on
    % ---------------------------------------------------------------------
    % Establish block boundary points -------------------------------------
    B = zeros(N+1,4);
    i = 1;
    while i <= N+1
        B(i,1) = R*cos(Bndry_blocks*(i-1));
        B(i,2) = R*sin(Bndry_blocks*(i-1));
        B(i,3) = RE*cos(Bndry_blocks*(i-1));
        B(i,4) = RE*sin(Bndry_blocks*(i-1));
        i = i+1;
    end
    % ---------------------------------------------------------------------
    % Set initial values and limits for the hinges and sliders ------------
    % Hinge H1
    H1 = 1;
    set(handles.slider1,'Min',0)
    set(handles.slider1,'Max',(N-1)/2-1)
    set(handles.slider1,'Value',H1-1)
    % Discretize the slider
    set(handles.slider1,'SliderStep',[(1/((N-1)/2-1)) (1/((N-1)/2-1))])
    % Hinge H2
    H2 = round((N+1)/4)-1; %this ensures and integer
    set(handles.slider2,'Min',H1)
    set(handles.slider2,'Max',(N-1)/2)
    set(handles.slider2,'Value',H2-1)
    % Discretize the slider
    set(handles.slider2,'SliderStep',[(1/((N-1)/2-H1)) (1/((N-1)/2-H1))])
    % Hinge H4
    H4 = N+1;
    set(handles.slider4,'Min',(N+1)/2+1)
    set(handles.slider4,'Max',N)
    set(handles.slider4,'Value',H4-1)
    % Discretize the slider
    set(handles.slider4,'SliderStep',...
        [(1/(H4-((N+1)/2))) (1/(H4-((N+1)/2)))])
    % Hinge H3
    H3 = round(3*(N+1)/4);
    set(handles.slider3,'Min',(N+1)/2)
    set(handles.slider3,'Max',H4-2)
    set(handles.slider3,'Value',H3-1)
    % Discretize the slider
    set(handles.slider3,'SliderStep',...
        [(1/(H4-((N+1)/2))) (1/(H4-((N+1)/2)))])
    % ---------------------------------------------------------------------
    % Plot the arch -------------------------------------------------------
    % Intrados
    plot(B(:,1),B(:,2),'k-')
    % Extrados
    plot(B(:,3),B(:,4),'k-')
    % Block boundaries
    for i = 1:N+1
        x = [B(i,1),B(i,3)];
        y = [B(i,2),B(i,4)];
        plot(x,y,'-','Color','k')
    end
    % set the axis
    axis manual
    YLim = get(gca,'YLim');
    YLim(2) = 1.25*YLim(2);
    YLim(1) = YLim(1)-0.01*(YLim(2)-YLim(1));
    set(gca,'YLim',YLim)
    % ---------------------------------------------------------------------
    % Hinge markers -------------------------------------------------------
    % Set hinge marker locations
    H_xy = [B(H1,3),B(H1,4);...
            B(H2,1),B(H2,2);...
            B(H3,3),B(H3,4);...
            B(H4,1),B(H4,2)];
    % Plot hinge markers
    hinge = ['1','2','3','4'];
    % Hinge 1
    handles.marker1 = plot(H_xy(1,1),H_xy(1,2),'or','MarkerSize',14,...
        'MarkerFaceColor','r');
    handles.marker1txt = text(H_xy(1,1),H_xy(1,2),hinge(1),...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    % Hinge 2
    handles.marker2 = plot(H_xy(2,1),H_xy(2,2),'or','MarkerSize',14,...
        'MarkerFaceColor','r');
    handles.marker2txt = text(H_xy(2,1),H_xy(2,2),hinge(2),...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    % Hinge 3
    handles.marker3 = plot(H_xy(3,1),H_xy(3,2),'or','MarkerSize',14,...
        'MarkerFaceColor','r');
    handles.marker3txt = text(H_xy(3,1),H_xy(3,2),hinge(3),...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    % Hinge 4
    handles.marker4 = plot(H_xy(4,1),H_xy(4,2),'or','MarkerSize',14,...
        'MarkerFaceColor','r');
    handles.marker4txt = text(H_xy(4,1),H_xy(4,2),hinge(4),...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    % ---------------------------------------------------------------------
    % Calculate center of mass (CM) location for full arch ----------------
    % CM angle
    CM_Aangle = pi/2;
    % CM radius
    CM_Arad = (4*sin((pi-0)/2)*(RE^3-R^3))/...
                            (3*(pi-0)*(RE^2-R^2)); %[m]
    % CM (x,y)
    CM_Axy(1,:) = [CM_Arad*cos(CM_Aangle),...
                            CM_Arad*sin(CM_Aangle)]; %[m,m]
    %Plot CM
    plot(CM_Axy(1,1),CM_Axy(1,2),'xk-')
    plot(CM_Axy(1,1),CM_Axy(1,2),'ok-')
    % ---------------------------------------------------------------------
    % Plot Load Label -----------------------------------------------------
    switch loadType
        % Point load
        case 1
            handles.PLt = text(B(H3,3),B(H3,4),{'\lambda_{P}';...
                    '\downarrow';' '},'HorizontalAlignment','center',...
                                            'VerticalAlignment','bottom');
            set(handles.PLt,'FontSize',18)
        % Horizontal acceleration
        case 2 
            handles.PLt = text(CM_Axy(1,1),CM_Axy(1,2),...
                '  \rightarrow \lambda_{a}*g','HorizontalAlignment',...
                                    'left','VerticalAlignment','middle');
            set(handles.PLt,'FontSize',18)   
    end
    % ---------------------------------------------------------------------
    % Thrust line plot check marker for handles update --------------------
    handles.TLPck = '0';
    % ---------------------------------------------------------------------
    % Run the Eval function -----------------------------------------------
    Eval(hObject, eventdata, handles)
    % ---------------------------------------------------------------------
end
% --------------------------- End RUN -------------------------------------
 
%% ----------------- EDIT TEXT & CHECKBOXES FOR INPUT ---------------------
% This section houses the create and callback functions for the input
% parameters. The input data is obtained from the RUN pushbutton so these
% parameters are unmodified. The one exception is the checkboxes for the 
% load types. There is a check in their callbacks to switch off the 
% checkbox of the other load case.
 
% N - number of blocks ----------------------------------------------------
function editText_N_Callback(hObject, eventdata, handles)
% hObject    handle to editText_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% --- Executes during object creation, after setting all properties.
function editText_N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editText_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
                                get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------

% R - clear span radius ---------------------------------------------------
function editText_R_Callback(hObject, eventdata, handles)
% hObject    handle to editText_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% --- Executes during object creation, after setting all properties.
function editText_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editText_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
                                get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------

% T_R - thickness to radius ratio -----------------------------------------
function editText_T_R_Callback(hObject, eventdata, handles)
% hObject    handle to editText_T_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
%- Executes during object creation, after setting all properties.
function editText_T_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editText_T_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
                                get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------

% Block depth -------------------------------------------------------------
function editText_Depth_Callback(hObject, eventdata, handles)
% hObject    handle to editText_Depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% --- Executes during object creation, after setting all properties.
function editText_Depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editText_Depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
                                get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------

% Block Density -----------------------------------------------------------
function editText_Density_Callback(hObject, eventdata, handles)
% hObject    handle to editText_Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% --- Executes during object creation, after setting all properties.
function editText_Density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editText_Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'),...
                                get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------

% PL checkbox -------------------------------------------------------------
% --- Executes on button press in checkbox_PL.
function checkbox_PL_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Check if H_acc checkbox is selected -------------------------------------
Hacc_chk = get(handles.checkbox_H_acc,'Value');
if Hacc_chk == 1
    set(handles.checkbox_H_acc,'Value',0);
end
% -------------------------------------------------------------------------
% H_acc checkbox ----------------------------------------------------------
% --- Executes on button press in checkbox_H_acc.
function checkbox_H_acc_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_H_acc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if PL checkbox is selected ----------------------------------------
PL_chk = get(handles.checkbox_PL,'Value');
if PL_chk == 1
    set(handles.checkbox_PL,'Value',0);
end
%--------------------------------------------------------------------------
% --------------- END EDIT TEXT & CHECKBOXES FOR INPUT --------------------

%% ------------------------- CLEAR BUTTON----------------------------------
% The clear button clears all data from the GUI.
% --- Executes on button press in clear_pushbutton.
function clear_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clear_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Clear the figure --------------------------------------------------------
cla(handles.figure_1, 'reset')
% -------------------------------------------------------------------------
% Clear all input parameters ----------------------------------------------
set(handles.editText_N, 'String',' '); 
set(handles.editText_R, 'String',' ');
set(handles.editText_T_R, 'String',' ');
set(handles.editText_Depth, 'String',' ');
set(handles.editText_Density, 'String',' ');
set(handles.checkbox_PL,'Value',0);
set(handles.checkbox_H_acc,'Value',0);
% -------------------------------------------------------------------------
% Clear the text on the hinge values --------------------------------------
set(handles.h1_text, 'String',' ');
set(handles.h2_text, 'String',' ');
set(handles.h3_text, 'String',' ');
set(handles.h4_text, 'String',' ');
% -------------------------------------------------------------------------
% Clear reaction values ---------------------------------------------------
set(handles.V1,'String',' ')
set(handles.V2,'String',' ')
set(handles.V3,'String',' ')
set(handles.V4,'String',' ')
set(handles.V5,'String',' ')
set(handles.V6,'String',' ')
set(handles.V7,'String',' ')
set(handles.V8,'String',' ')
% -------------------------------------------------------------------------
% Clear lambda value ------------------------------------------------------
set(handles.lamda, 'String',' ');
% -------------------------------------------------------------------------
% Set sliders to the minimum value ----------------------------------------
reset_slider1 = get(handles.slider1,'Min');
set(handles.slider1, 'Value', reset_slider1);
reset_slider2 = get(handles.slider2,'Min');
set(handles.slider2, 'Value', reset_slider2);
reset_slider3 = get(handles.slider3,'Min');
set(handles.slider3, 'Value', reset_slider3);
reset_slider4 = get(handles.slider4,'Min');
set(handles.slider4, 'Value', reset_slider4);
% -------------------------------------------------------------------------
% Update all handles and clear all ----------------------------------------
guidata(hObject, handles);
clear all   
% --------------------------- END CLEAR -----------------------------------
 
%% ----------------------EVALUATION FUNCTION ------------------------------
% The Eval function is where all of the analysis calculations are
% performed. This function is called by the RUN button and by each slider
% callback. When it is called it calculates the element properties of the
% three segments of the arch that form the mechanism. Each elements CM
% location, weight, and lever arms are established and supplied to the
% matrix form of the equilibrium equations. These equations are developed 
% on the assumption of rigid elements between perfect hinges and solved to
% determine the collapse load and reactions at the four hinges. The
% resulting values are then chacked to ensure that the collapse muliplier
% is positive and the no tension assumption is upheld. It this is true the
% results are displayed. Finally, using the reactions at H1 the thrust line
% is determined by calculating the eccentricity from the centerline of the
% arch required to maintain zero moment at H1. This line is then plotted.
%
function Eval(hObject, eventdata, handles)
% hObject    handle to run_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Obtain data from handles ------------------------------------------------
% N - Number of blocks
N = str2num(get(handles.editText_N, 'String'));
% R - Internal radius
R = str2num(get(handles.editText_R,'String')); % [m]
% T_R - Thickness to radius ratio
T_R = str2num(get(handles.editText_T_R,'String')); 
% Block depth
Depth_blocks=str2num(get(handles.editText_Depth,'String')); %[m]
% Block density
Den_blocks=str2num(get(handles.editText_Density,'String')); %[kg/m^3]
% Point load checkbox value
Pl=get(handles.checkbox_PL,'Value');
% Horizontal acceleration checkbox value
H_acc=get(handles.checkbox_H_acc,'Value');
% Hinge 1
H1=get(handles.slider1,'Value')+1;
% Hinge 2
H2=get(handles.slider2,'Value')+1;
% Hinge 3
H3=get(handles.slider3,'Value')+1;
% Hinge 4
H4=get(handles.slider4,'Value')+1;
% -------------------------------------------------------------------------
% Base calculations -------------------------------------------------------
% Thickness
T = T_R*R; %[m]
% Extrados Radius
RE = R+T; %[m]
% Bndry_blocks - The step between block boundaries
Bndry_blocks = pi/N;
% -------------------------------------------------------------------------
% Set load type -----------------------------------------------------------
% 1 for point load
if Pl == 1
    loadType = 1;
% 2 for horizontal acceleration
elseif H_acc == 1
    loadType = 2;
end
% -------------------------------------------------------------------------
% Gravity -----------------------------------------------------------------
g = 9.81; %[m/s^2]
% -------------------------------------------------------------------------
% Establish block boundary points -----------------------------------------
B = zeros(N+1,4);
i = 1;
while i <= N+1
    B(i,1) = R*cos(Bndry_blocks*(i-1));
    B(i,2) = R*sin(Bndry_blocks*(i-1));
    B(i,3) = RE*cos(Bndry_blocks*(i-1));
    B(i,4) = RE*sin(Bndry_blocks*(i-1));
    i = i+1;
end
% -------------------------------------------------------------------------
% Set single hinge variable for loop calculations -------------------------
H_all = [H1,H2,H3,H4];
H_theta = Bndry_blocks*(H_all-1);
% -------------------------------------------------------------------------
% Designate joint location for point load ---------------------------------
if loadType == 1
    Pl_Theta = Bndry_blocks*H3;
    PlRad = RE;
end
% -------------------------------------------------------------------------
% Establish the total Arch CM location ------------------------------------
CM_Aangle = pi/2;
CM_Arad = (4*sin((pi-0)/2)*(RE^3-R^3))/...
                    (3*(pi-0)*(RE^2-R^2)); %[m]
CM_Axy(1,:) = [CM_Arad*cos(CM_Aangle),...
                        CM_Arad*sin(CM_Aangle)]; %[m,m]
% -------------------------------------------------------------------------
% Calculate each elements CM, area, and gravitational force ---------------
% (assuming uniform density and dimensions)
% Set variables
CM_angle=zeros(1,3);
CM_rad=zeros(1,3);
Area=zeros(1,3);
Fg=zeros(1,3);
CM_xy=zeros(3,2);
for i=1:3
    % Establish center of mass coordinates
    CM_angle(i) = (H_theta(i+1)+H_theta(i))/2;
    CM_rad(i) = (4*sin((H_theta(i+1)-H_theta(i))/2)*(RE^3-R^3))/...
                (3*(H_theta(i+1)-H_theta(i))*(RE^2-R^2)); %[m]
    CM_xy(i,:) = [CM_rad(i)*cos(CM_angle(i)),...
                    CM_rad(i)*sin(CM_angle(i))]; %[m,m]
    % Establish Elements Cross sectional area
    Area(i) = ((H_theta(i+1)-H_theta(i))/2)*(RE^2-R^2); %[m^2]
    % Calculate gravitational force 
    Fg(i) = (Den_blocks*g*Depth_blocks*Area(i))/1000; % kN      
end
% -------------------------------------------------------------------------
% Establish Cartesian Coordinates for hinges ------------------------------
H_xy = [B(H1,3),B(H1,4);...
        B(H2,1),B(H2,2);...
        B(H3,3),B(H3,4);...
        B(H4,1),B(H4,2)];
% -------------------------------------------------------------------------
% Display (x,y) coordinates for the hinge locations -----------------------
% Limit numbers to 4 digits for display
HxyC = chop(H_xy,4);
% Hinge 1
H1xy = ['(',num2str(HxyC(1,1)),',',num2str(HxyC(1,2)),')']; %[m]
set(handles.h1_text,'String',H1xy)
% Hinge 2
H2xy = ['(',num2str(HxyC(2,1)),',',num2str(HxyC(2,2)),')']; %[m]
set(handles.h2_text,'String',H2xy)
% Hinge 3
H3xy = ['(',num2str(HxyC(3,1)),',',num2str(HxyC(3,2)),')']; %[m]
set(handles.h3_text,'String',H3xy)
% Hinge 4
H4xy = ['(',num2str(HxyC(4,1)),',',num2str(HxyC(4,2)),')']; %[m]
set(handles.h4_text,'String',H4xy)
% -------------------------------------------------------------------------
%% --------------------- Equilibrium Problem ------------------------------
% Balance Equation [BC][V]=[Q]
% [BC] = Balance equation variable constants
% [V] = Variables [h1;v1;h2;v2;h3;v3;h4;v4;Lambda]
% h - horizontal reaction
% v - vertical reaction
% [Q] = Constants (gravity's contribution)
% Calculate moment arms for equilibrium equation --------------------------
% Moment arms are incorporated into [BC]
switch loadType
    case 1 % Point Load
    % Row1 = ME1@h1[H2,V2,Lambda,gravity]
        M = [(H_xy(2,2)-H_xy(1,2)),...
            (H_xy(1,1)-H_xy(2,1)),...
             0,...
            (Fg(1)*(H_xy(1,1)-CM_xy(1,1)));
    % Row2 = ME2@h2[H3,V3,Lambda,gravity]
            (H_xy(3,2)-H_xy(2,2)),...
            (H_xy(2,1)-H_xy(3,1)),...
            (H_xy(2,1)-H_xy(3,1)),...
            (Fg(2)*(H_xy(2,1)-CM_xy(2,1)));
    % Row3 = ME3@h3[H4,V4,Lambda,gravity]
            (H_xy(3,2)-H_xy(4,2)),...
            (H_xy(3,1)-H_xy(4,1)),...
            0,...
            (Fg(3)*(H_xy(3,1)-CM_xy(3,1)))];
    case 2  % Horizontal Acceleration
    % Row1 = ME1@h1[H2,V2,Lambda,gravity]
        M = [(H_xy(2,2)-H_xy(1,2)),...
            (H_xy(1,1)-H_xy(2,1)),...
            (Fg(1)*(CM_xy(1,2)-H_xy(1,2))),...
            (Fg(1)*(H_xy(1,1)-CM_xy(1,1)));
    % Row2 = ME2@h2[H3,V3,Lambda,gravity]
            (H_xy(3,2)-H_xy(2,2)),...
            (H_xy(2,1)-H_xy(3,1)),...
            (Fg(2)*(CM_xy(2,2)-H_xy(2,2))),...
            (Fg(2)*(H_xy(2,1)-CM_xy(2,1)));
    % Row3 = ME3@h3[H4,V4,Lambda,gravity]
            (H_xy(3,2)-H_xy(4,2)),...
            (H_xy(3,1)-H_xy(4,1)),...
            (Fg(3)*(H_xy(3,2)-CM_xy(3,2))),...
            (Fg(3)*(H_xy(3,1)-CM_xy(3,1)))];
end
% -------------------------------------------------------------------------
% Establish [BC] and [Q] --------------------------------------------------
switch loadType
    case 1
        % Row1:Row3 = Fx@E1;Fy@E1;M@E1
        BC = [-1,0,1,0,0,0,0,0,0;...
            0,1,0,-1,0,0,0,0,0;...
            0,0,-M(1,1),M(1,2),0,0,0,0,0;...
        % Row4:Row6 = Fx@E2;Fy@E2;M@E2
            0,0,-1,0,1,0,0,0,0;...
            0,0,0,1,0,1,0,0,-1;...
            0,0,0,0,M(2,1),M(2,2),0,0,-M(2,3);...
        % Row7:Row9 = Fx@E3;Fy@E3;M@E3
            0,0,0,0,-1,0,1,0,0;...
            0,0,0,0,0,-1,0,1,0;...
            0,0,0,0,0,0,M(3,1),-M(3,2),0];
        % end [BC]
        Q = [0;Fg(1);-M(1,4);0;Fg(2);M(2,4);0;Fg(3);-M(3,4)];
    case 2
        % Row1:Row3 = Fx@E1;Fy@E1;M@E1
        BC = [-1,0,1,0,0,0,0,0,Fg(1);...
            0,1,0,-1,0,0,0,0,0;...
            0,0,-M(1,1),M(1,2),0,0,0,0,-M(1,3);...
        % Row4:Row6 = Fx@E2;Fy@E2;M@E2
            0,0,-1,0,1,0,0,0,Fg(2);...
            0,0,0,1,0,1,0,0,0;...
            0,0,0,0,M(2,1),M(2,2),0,0,M(2,3);...
        % Row7:Row9 = Fx@E3;Fy@E3;M@E3
            0,0,0,0,-1,0,1,0,Fg(3);...
            0,0,0,0,0,-1,0,1,0;...
            0,0,0,0,0,0,M(3,1),-M(3,2),M(3,3)];
        % end BC
        Q = [0;Fg(1);-M(1,4);0;Fg(2);M(2,4);0;Fg(3);-M(3,4)];
end
% -------------------------------------------------------------------------
% Calculate collapse load and reactions -----------------------------------
% Use the reverse division procedure in MATLAB for calculation
V = BC\Q;
% ------------------------ End Equilibrium Problem ------------------------
 
%% ------------------- THRUST LINE CALCULATION ----------------------------
% This section calculates the coefficients of the thrust line equation.
% The thrust line is a load path that follows the concentration of
% compressive forces. From the determined reactions V(1:8) and the
% collapse load V(9) the thrust line location is determined by
% calculating the eccentricity required to maintain a zero moment at
% hinge H1.
switch loadType
    case 1
        % Point load calculation is divided into two segments to account
        % for the inclusion of the point load.
        % Set thrust line data set 1 and increment (H1 to H3) -------------
        Tl1 = zeros(200,2);
        Theta_delta = (H_theta(3)-H_theta(1))/200;
        % -----------------------------------------------------------------
        for i = 1:200
            % Calculate sub CM --------------------------------------------
            % (identified with CMhat)
            CMhat_angle = H_theta(1) + Theta_delta*(i)/2;
            CMhat_rad = (4*sin((Theta_delta*(i))/2)*(RE^3-R^3))/...
                        (3*(Theta_delta*(i))*(RE^2-R^2)); %[m]
            CMhat_xy = [CMhat_rad*cos(CMhat_angle),...
                        CMhat_rad*sin(CMhat_angle)]; %[m,m]
            % -------------------------------------------------------------
            % Calculate cross sectional area ------------------------------
            Ahat = ((Theta_delta*(i))/2)*(RE^2-R^2); %[m^2]
            % -------------------------------------------------------------
            % Calculate gravitational force -------------------------------
            Fghat = (Den_blocks*g*Depth_blocks*Ahat)/1000; %[kN]
            % -------------------------------------------------------------
            % Reactions ---------------------------------------------------
            % Horizontal
            H_hat = V(1); %[kN]
            % Vertical
            V_hat = V(2) - Fghat; %[kN]
            % -------------------------------------------------------------
            % Establish moment lever arms for the moment taken about H1 ---
            DxCM_hat = H_xy(1,1)- CMhat_xy(1); %[m]
            DyCM_hat = CMhat_xy(2) - H_xy(1,2); %[m]
            % -------------------------------------------------------------
            % Calculate eccentricity from equilibrium equations -----------
            % e_R = (PI_R-(A_R+B_R)*Rcl)/(A_R+B_R)
            PI_R = Fghat*((DxCM_hat))+V_hat*H_xy(1,1)+H_hat*H_xy(1,2);
            A_R = H_hat*sin(H_theta(1) + Theta_delta*(i));
            B_R = V_hat*cos(H_theta(1) + Theta_delta*(i));
            e_R = (PI_R-(A_R+B_R)*(R+T/2))/(A_R+B_R); %[m]
            % -------------------------------------------------------------
            % Set cartesian coordinates -----------------------------------
            x_hat = (R+T/2+e_R)*cos(H_theta(1) + Theta_delta*(i)); %[m]
            y_hat = (R+T/2+e_R)*sin(H_theta(1) + Theta_delta*(i)); %[m]
            % -------------------------------------------------------------
            % Set thrust line data point ----------------------------------
            Tl1(i,:) = [x_hat,y_hat]; %[m,m]
            % -------------------------------------------------------------
        end
        % Set thrust line data set 2 and increment (H3 to H4) -------------
        Tl2 = zeros(200,2);
        Theta_delta = (H_theta(4)-H_theta(3))/200;
        % -----------------------------------------------------------------
        for i = 1:200
            % Calculate sub CM --------------------------------------------
            % (identified with CMhat)
            CMhat_angle = (H_theta(3)+Theta_delta*(i)+H_theta(1))/2;
            CMhat_rad = (4*sin(((H_theta(3)+Theta_delta*(i))...
                               -H_theta(1))/2)*(RE^3-R^3))/(3*...
                                    ((H_theta(3)+Theta_delta*(i))...
                                            - H_theta(1))*(RE^2-R^2)); %[m]
            CMhat_xy = [CMhat_rad*cos(CMhat_angle),...
                        CMhat_rad*sin(CMhat_angle)]; %[m,m]
            % -------------------------------------------------------------
            % Calculate cross sectional area ------------------------------
            Ahat = (((H_theta(3)+Theta_delta*(i))...
                                         -H_theta(1))/2)*(RE^2-R^2); %[m^2]
            % -------------------------------------------------------------
            % Calculate gravitational force -------------------------------
            Fghat = (Den_blocks*g*Depth_blocks*Ahat)/1000; %[kN]
            % -------------------------------------------------------------
            % Reactions ---------------------------------------------------
            % Horizontal
            H_hat = V(1); %[kN]
            % Vertical
            V_hat = V(2)- V(9) - Fghat; %[kN]
            % -------------------------------------------------------------
            % Establish moment lever arms for the moment taken about H1 ---
            DxCM_hat = H_xy(1,1)- CMhat_xy(1); %[m]
            DyCM_hat = CMhat_xy(2) - H_xy(1,2); %[m]
            % -------------------------------------------------------------
            % Calculate eccentricity from equilibrium equations -----------
            % e_R = (PI_R-(A_R+B_R)*Rcl)/(A_R+B_R)
            PI_R = Fghat*((DxCM_hat))+V_hat*H_xy(1,1)+H_hat*H_xy(1,2)...
                                               +V(9)*(H_xy(1,1)-H_xy(3,1));
            A_R = H_hat*sin(Theta_delta*(i)+H_theta(3));
            B_R = V_hat*cos(Theta_delta*(i)+H_theta(3));
            e_R = (PI_R-(A_R+B_R)*(R+T/2))/(A_R+B_R); %[m]
            % -------------------------------------------------------------
            % Set cartesian coordinates -----------------------------------
            x_hat = (R+T/2+e_R)*cos(Theta_delta*(i)+H_theta(3)); %[m]
            y_hat = (R+T/2+e_R)*sin(Theta_delta*(i)+H_theta(3)); %[m]
            % -------------------------------------------------------------
            % Set thrust line data point ----------------------------------
            Tl2(i,:) = [x_hat,y_hat]; %[m,m]
            % -------------------------------------------------------------
        end 
    case 2 
        % Horizontal acceleration -----------------------------------------
        % Horizontal acceleration is applied smoothly, only 1 thrust line
        % required
        % Set thrust line data set 1 and increment (H1 to H4) -------------
        Tl1 = zeros(200,2);
        Theta_delta = (H_theta(4)-H_theta(1))/200;
        % -----------------------------------------------------------------
        for i = 1:200
            % Calculate sub CM --------------------------------------------
            % (identified with CMhat)
            CMhat_angle = H_theta(1) + Theta_delta*(i)/2;
            CMhat_rad = (4*sin((Theta_delta*(i))/2)*(RE^3-R^3))/...
                        (3*(Theta_delta*(i))*(RE^2-R^2)); %[m]
            CMhat_xy = [CMhat_rad*cos(CMhat_angle),...
                        CMhat_rad*sin(CMhat_angle)]; %[m,m]
            % -------------------------------------------------------------
            % Calculate cross sectional area ------------------------------
            Ahat = ((Theta_delta*(i))/2)*(RE^2-R^2); %[m^2]
            % -------------------------------------------------------------
            % Calculate gravitational force -------------------------------
            Fghat = (Den_blocks*g*Depth_blocks*Ahat)/1000; %[kN]
            % -------------------------------------------------------------
            % Reactions ---------------------------------------------------
            % Horizontal
            H_hat = V(1) - Fghat*V(9); %[kN]
            % Vertical
            V_hat = V(2) - Fghat; %[kN]
            % -------------------------------------------------------------
            % Establish moment lever arms for the moment taken about H1 ---
            DxCM_hat = H_xy(1,1)- CMhat_xy(1); %[m]
            DyCM_hat = CMhat_xy(2) - H_xy(1,2); %[m]
            % -------------------------------------------------------------
            % Calculate eccentricity from equilibrium equations -----------
            % e_R = (PI_R-(A_R+B_R)*Rcl)/(A_R+B_R)
            PI_R = Fghat*((DxCM_hat-V(9)*DyCM_hat))...
                                        +V_hat*H_xy(1,1)+H_hat*H_xy(1,2);
            A_R = H_hat*sin(H_theta(1) + Theta_delta*(i));
            B_R = V_hat*cos(H_theta(1) + Theta_delta*(i));
            e_R = (PI_R-(A_R+B_R)*(R+T/2))/(A_R+B_R); %[m]
            % -------------------------------------------------------------
            % Set cartesian coordinates -----------------------------------
            x_hat = (R+T/2+e_R)*cos(H_theta(1) + Theta_delta*(i)); %[m]
            y_hat = (R+T/2+e_R)*sin(H_theta(1) + Theta_delta*(i)); %[m]
            % -------------------------------------------------------------
            % Set thrust line data point ----------------------------------
            Tl1(i,:) = [x_hat,y_hat]; %[m]
            % -------------------------------------------------------------
        end
end
% --------------------- End Thrust Line Calculation -----------------------
 
%% ----------------- Collapse Load and Reaction Display -------------------
% Now that the calculations are all performed, the data needs to be checked
% for admissibility. The two criteria are if the collapse load is negative
% and if the reactions at the hinges violate the no tension rule. The no
% tension rule is examined by comparing the angle of the net reaction
% against the boundary surface.
% Setting check parameter for thrust line plots ---------------------------
TLPck = str2num(handles.TLPck);
% -------------------------------------------------------------------------
% Setup check for no-tension violation of the normal force ----------------
% hinge angles
thetaH = H_theta;
% Tension check parameter
negChk = 0;
% Set the reaction angle parameter 
thetaN = zeros(1,4);
% Set default case check value
caseChk = 0;
% -------------------------------------------------------------------------
for i = 2:2:8
    % Case check ----------------------------------------------------------
    % 1 - Negative vertical reaction
    if V(i) < 0
        caseChk = 1;
    % 2- Negative horizontal reaction
    elseif V(i-1) < 0
        caseChk = 2;
    % 3 - Negative vertical and horizontal reactions
    elseif V(i) < 0 && V(i-1) < 0
        caseChk = 3;
    end
    % ---------------------------------------------------------------------
    % Set reaction angle values based on case type ------------------------
    switch caseChk
        case 1
            thetaN(i/2) = atan((-V(i))/V(i-1));
        case 2
            thetaN(i/2) = atan(V(i)/(-V(i-1)));
        case 3
            thetaN(i/2) = 5*pi;
    end
    % ---------------------------------------------------------------------
    % Perform tension checks ----------------------------------------------
    % Horizontal block joint
    if thetaH(i/2) == 0 || thetaH(i/2) == pi
        if caseChk == 1
            negChK = 1;
        end
    % First quadrant
    elseif thetaH(i/2) < (pi/2)
        switch caseChk
            case 1
                if thetaN(i/2) >= thetaH(i/2)
                    negChk = 1;
                end
            case 2
                if thetaN(i/2) <= thetaH(i/2)
                    negChk = 1;
                end
            case 3
                negChk = 1;
        end
    % Second quadrant
    elseif thetaH(i/2) > (pi/2)
        switch caseChk
            case 1
                if thetaN(i/2) >= (pi-thetaH(i/2))
                    negChk = 1;
                end
            case 2
                if thetaN(i/2) <= (pi-thetaH(i/2))
                    negChk = 1;
                end
            case 3
                negChk = 1;
        end
    end
end
% NOTE: Vertical joint is not checked due to the odd number of blocks
% requirement.
% -------------------------------------------------------------------------
% Set collapse load to 'Not admissible' and reactions blank for violations
if V(9) < 0 || negChk == 1
    Tot = ['Not admissible'];
    % Display reaction values
    V1 = [];
    V2 = [];
    V3 = [];
    V4 = [];
    V5 = [];
    V6 = [];
    V7 = [];
    V8 = [];
    % Clear thrust line data if it exists
    if TLPck == 1 && loadType == 1
        set(handles.TLptR,'XData',[],'YData',[])
        set(handles.TLptL,'XData',[],'YData',[])
    elseif TLPck == 1 && loadType == 2
        set(handles.TLptR,'XData',[],'YData',[])
    end
else
    % Set collapse load and reactions for display to four digits ----------
    % Collapse load
    tot=chop(V(9),4);
    % Reactions and load in string format
    V1 = num2str(chop(V(1),4)); %[kN]
    V2 = num2str(chop(V(2),4)); %[kN]
    V3 = num2str(chop(V(3),4)); %[kN]
    V4 = num2str(chop(V(4),4)); %[kN]
    V5 = num2str(chop(V(5),4)); %[kN]
    V6 = num2str(chop(V(6),4)); %[kN]
    V7 = num2str(chop(V(7),4)); %[kN]
    V8 = num2str(chop(V(8),4)); %[kN]
    tot=num2str(tot);
    % Setup collapse load display based on loading type
    switch loadType
         case 1 % Point Load
             Tot = [tot,' kN'];
         case 2 % Horizontal Acceleration 
            Tot = tot;
    end
end
% -------------------------------------------------------------------------
% Display collapse load and reactions -------------------------------------
set(handles.lamda,'String',Tot);
set(handles.V1,'String',V1);
set(handles.V2,'String',V2);
set(handles.V3,'String',V3);
set(handles.V4,'String',V4);
set(handles.V5,'String',V5);
set(handles.V6,'String',V6);
set(handles.V7,'String',V7);
set(handles.V8,'String',V8);
% -------------------------------------------------------------------------
% ---------------- End Collapse Load and Reaction Display -----------------
%% ------------------------ Plot Thrust Line ------------------------------
% Plot the trust line if admissible solution ------------------------------
if V(9) >= 0 && negChk == 0
    switch loadType
        % Point Load
        case 1
           % If the thrust line plot handles do not exist -----------------
           if TLPck == 0
               % Set figure axis as current
                axes(handles.figure_1); 
                hold on
                axis manual
                % Plot left and right thrust lines
                handles.TLptR = plot(Tl1(:,1),Tl1(:,2),'b');
                handles.TLptL = plot(Tl2(:,1),Tl2(:,2),'b');
                handles.TLPck = '1';
           % -------------------------------------------------------------
           % If the thrust line plot handles exist -----------------------
           else
                set(handles.TLptR,'XData',Tl1(:,1),'YData',Tl1(:,2))
                set(handles.TLptL,'XData',Tl2(:,1),'YData',Tl2(:,2))   
           end 
           % --------------------------------------------------------------
        % Horizontal Acceleration
        case 2 
            % If the thrust line plot handles do not exist ----------------
            if TLPck == 0
                % Set figure axis as current
                axes(handles.figure_1); 
                hold on
                axis manual
                % Plot thrust line
                handles.TLptR = plot(Tl1(:,1),Tl1(:,2),'b');
                set(handles.TLptR,'Tag','test')
                handles.TLPck = '1';
            % -------------------------------------------------------------
            % If the thrust line plot handles exist -----------------------
            else
                set(handles.TLptR,'XData',Tl1(:,1),'YData',Tl1(:,2))
            end
            % -------------------------------------------------------------
        hold off
    end
end  
% Update the handles
guidata(hObject, handles);
% ---------------------- End Plot Thrust Line -----------------------------
% ----------------------END EVALUATION FUNCTION ---------------------------
% ---------------------END KCLCalculator.m -----------------------


