function varargout = GUI_Filter(varargin)
% GUI_FILTER MATLAB code for GUI_Filter.fig
%      GUI_FILTER, by itself, creates a new GUI_FILTER or raises the existing
%      singleton*.
%
%      H = GUI_FILTER returns the handle to a new GUI_FILTER or the handle to
%      the existing singleton*.
%
%      GUI_FILTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FILTER.M with the given input arguments.
%
%      GUI_FILTER('Property','Value',...) creates a new GUI_FILTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Filter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Filter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Filter

% Last Modified by GUIDE v2.5 10-Nov-2011 18:40:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Filter_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Filter_OutputFcn, ...
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
end

% --- Executes just before GUI_Filter is made visible.
function GUI_Filter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Filter (see VARARGIN)

% Choose default command line output for GUI_Filter
handles.output = hObject;
% Give the acces to the figure Toolbar
set(hObject,'toolbar','figure');

handles.Hnum = 0; 
handles.Hden = 0;
handles.pole = 0;
handles.zero = 0;
handles.gain = 0;
handles.z = 0;
handles.x = 0;
handles.y = 0;
handles.dtftLinIndex = 0;
handles.dtftIndex = 0;
handles.notDTFTIndex = 0;
handles.R = 0; 
handles.H = 0;

handles.w_c = 0;
handles.w_c2 = 0;
handles.wprime_c = 0;
handles.wprime_c2 = 0;

handles.w = 0;
handles.wTransform =0; 

handles.type = 0; 
handles.alpha = 0;
handles.beta = 0;
handles.ZtransNum = 0;
handles.ZtransDen = 0; 


handles.HnumTrans = 0;
handles.HdenTrans = 0;

handles.dtftLinIndexTrans = 0;
handles.dtftIndexTrans = 0;
handles.notDTFTIndexTrans = 0;
handles.RTrans = 0; 
handles.HTrans = 0;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Filter wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Filter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transfer  Function%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in validate_tf_pushbutton.
function [handles]= validate_tf_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to validate_tf_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axis = get(handles.root_locus_axis);
hold off; 
[handles] = getH(hObject,handles);

set(handles.w_c_edittext,'String', num2str(handles.wco));
guidata(hObject, handles);
zplane(handles.Hnum,handles.Hden); 
hold on; 
set(handles.freqInfo_statictext, 'String', '\alpha =  and  \Beta = ');
set(handles.num_freq_tran_statictext, 'String', '0');
set(handles.den_freq_tran_statictext, 'String', '0'); 
    
end


function num_edittext_Callback(hObject, eventdata, handles)
% hObject    handle to num_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_edittext as text
%        str2double(get(hObject,'String')) returns contents of num_edittext as a double


%num = str2num(get(handles.num_edittext,'String')); 
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function num_edittext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles);
end

function den_edittext_Callback(hObject, eventdata, handles)
% hObject    handle to den_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of den_edittext as text
%        str2double(get(hObject,'String')) returns contents of den_edittext as a double
%den = str2num(get(handles.den_edittext,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function den_edittext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to den_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in num_den_radiobutton.
function num_den_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to num_den_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of num_den_radiobutton

% Allow only the 'pole/zero' mode or the 'num/den' mode
if (get(handles.pole_zero_radiobutton,'Value'))
    set(handles.pole_zero_radiobutton,'Value',0); 
end
   
guidata(hObject, handles);
end

% --- Executes on button press in pole_zero_radiobutton.
function pole_zero_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to pole_zero_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pole_zero_radiobutton

% Allow only the 'pole/zero' mode or the 'num/den' mode
if (get(handles.num_den_radiobutton,'Value'))
    set(handles.num_den_radiobutton,'Value',0); 
end
   
guidata(hObject, handles);
end

function pole_edittext_Callback(hObject, eventdata, handles)
% hObject    handle to pole_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pole_edittext as text
%        str2double(get(hObject,'String')) returns contents of pole_edittext as a double

%pole = str2num(get(handles.pole_edittext, 'String'));

if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function pole_edittext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pole_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function zero_edittext_Callback(hObject, eventdata, handles)
% hObject    handle to zero_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zero_edittext as text
%        str2double(get(hObject,'String')) returns contents of zero_edittext as a double

%zero = str2num(get(handles.zero_edittext, 'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function zero_edittext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zero_edittext (see GCBO)pole_zero
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function gain_edittext_Callback(hObject, eventdata, handles)
% hObject    handle to gain_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gain_edittext as text
%        str2double(get(hObject,'String')) returns contents of gain_edittext as a double
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function gain_edittext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gain_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function [num den] = zpg2nd(zero, pole, gain)
    num = gain*poly(zero);
    den = poly(pole); 
end

function  [zero pole gain ] = nd2zpg(num,den)
    den2 = den./den(1); 
    num2 = num./den(1); 
    gain = num2(1);
    zero = roots(num2)'; 
    pole = roots(den2)'; 
end  


% --- Executes on button press in freq_tran_pushbutton.
function [handles]=freq_tran_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to freq_tran_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
alpha = 0;
beta = 1; 

z = handles.z; 
ZtransNum = 0;
ZtransDen = 0; 

Hnum = handles.Hnum;
Hden = handles.Hden; 

    %Recuparation of the different cut-off frequencies
    w_c = str2num(get(handles.w_c_edittext, 'String'));
    wprime_c  = str2num(get(handles.wprime_c_edittext, 'String'));
    % And only in case of Band Pass or Stop Band 
    wprime_c2 = str2num(get(handles.wprime_c2_edittext, 'String'));
    
    
    type = get(handles.transform_type_popupmenu,'Value');
    numOrder = length(Hnum) - 1;
    denOrder = length(Hden) - 1;  
   
    % In case we apply two frequency transformations in a row, we want to
    % keep the original zplane.
    
    axis = get(handles.root_locus_axis);
    hold off;
    zplane(Hnum,Hden);
    hold on;   
    
    % To get HTrans
    
    ii = 0;
    
   % hold on; 
    for fToTransformTo = linspace(w_c,wprime_c,1)
        ii = ii+1;
            %%%
        switch type

            case 2 % Low Pass

                alpha = sin((w_c-fToTransformTo)/2)/sin((w_c+fToTransformTo)/2);
                beta = 'betaLP';

                zTransformedNum = [1 -alpha];
                zTransformedDen = [-alpha 1]; 



            case 3 % High Pass

                alpha = -cos((w_c+fToTransformTo)/2)/cos((w_c-fToTransformTo)/2);
                beta = 'betaHP';

                zTransformedNum = [-1 -alpha];
                zTransformedDen = [alpha 1];
           

            case 4 % Band Pass

                alpha = cos((wprime_c2+fToTransformTo)/2)/cos((wprime_c2-fToTransformTo)/2);
                beta = cot((wprime_c2-fToTransformTo)/2)*tan(w_c/2); 

                zTransformedNum = [-1,(2*alpha*beta)/(beta+1), -(beta-1)/(beta+1)];
                zTransformedDen = [(beta-1)/(beta+1), (-2*alpha*beta)/(beta+1) 1];


            case 5 % Stop Band

                alpha =  cos((wprime_c2+fToTransformTo)/2)/cos((fToTransformTo-wprime_c2)/2);
                beta = tan((wprime_c2-fToTransformTo)/2)*tan(w_c/2); 

                zTransformedNum = [1,-(2*alpha*beta)/(beta+1), (beta-1)/(beta+1)];
                zTransformedDen = [-(beta-1)/(beta+1), (-2*alpha)/(beta+1) 1];


            otherwise 
                alpha = 0;
                beta = 1; 
                numOrder = length(Hnum) - 1;
                denOrder = length(Hden) - 1;
        end
        
        % Let's show the value of the parameters alpha and beta
        set(handles.freqInfo_statictext, 'String', ['alpha = ' num2str(alpha) ' and  beta = ' num2str(beta)]);
    
        zTransformedNumOrder = length(zTransformedNum) - 1;
        zTransformedDenOrder = length(zTransformedDen) - 1;
        zTransNumSum = 0;
        zTransDenSum = 0;
        
        wTransformNum = 0;
        wTransformDen = 0;
        w = linspace(-pi,pi,500);
        for nOrder = 0 : zTransformedNumOrder
            zTransNumSum = zTransNumSum + zTransformedNum(nOrder+1)*z.^(zTransformedNumOrder - nOrder); 
            wTransformNum = wTransformNum + zTransformedNum(nOrder+1)*w.^(zTransformedNumOrder - nOrder);
        end
        for dOrder = 0 : zTransformedDenOrder
            zTransDenSum = zTransDenSum + zTransformedDen(dOrder+1)*z.^(zTransformedDenOrder - dOrder); 
            wTransformDen = wTransformDen + zTransformedDen(dOrder+1)*w.^(zTransformedDenOrder - dOrder);
        end
        w_1 = w;
         zTrans = zTransNumSum ./ zTransDenSum;
         wTransform = wTransformNum./wTransformDen;
        % wTransform(find((wTransform)>50))=50;
        % wTransform(find((wTransform)<-50))=-50;

         
        % Calculate H
        numSum = 0;
        denSum = 0;
        for numPower = 0 : numOrder
            numSum = numSum + Hnum(numPower+1)*zTrans.^(numOrder - numPower); 
        end
        for denPower = 0 : denOrder
            denSum = denSum + Hden(denPower+1)*zTrans.^(denOrder - denPower);
        end

        HTrans = numSum./denSum; 
        %zeroIndeces = find(HTrans == 0);

        % Find coefficients of z that are at the unit circle
        dtftLinIndexTrans = find((abs(z) > 0.999) & (abs(z) < 1.001));
        dtftIndexTrans = ind2sub(size(z),dtftLinIndexTrans);
        dtftTrans = HTrans(dtftIndexTrans);
        notDTFTIndexTrans = find(abs(z) > 1.025 | abs(z) < .975);
        w = angle(z(dtftIndexTrans));
        posFIndexTrans = find(w>0);
        
        % Find the numerator and denumerator of the new transfer function 
        % Band Pass and Stop Band : the degree is one higher
        %
        if (type > 3)
            [HnumTrans,HdenTrans] = invfreqz(dtftTrans(posFIndexTrans),w(posFIndexTrans),numOrder+1,denOrder+1);
        else
            [HnumTrans,HdenTrans] = invfreqz(dtftTrans(posFIndexTrans),w(posFIndexTrans),numOrder,denOrder);
        end

        % determine the zeros and poles with the numerators and denominators
        [zeroTrans, poleTrans, gainTrans] = nd2zpg(HnumTrans,HdenTrans);
 
        %figure(1); 
        plot(real(zeroTrans),imag(zeroTrans),'o','Color',[.1*ii 0 0]);hold on; 
        plot(real(poleTrans),imag(poleTrans),'x','Color',[.1*ii 0 0]);hold on; 
    end
    %hold off; 
    % Show the new coefficiant (2 decimals)
    set(handles.num_freq_tran_statictext, 'String', num2str((round(HnumTrans.*100))/100));
    set(handles.den_freq_tran_statictext, 'String', num2str((round(HdenTrans.*100))/100)); 
    
    RTrans = zeros(size(z));
    RTrans(dtftIndexTrans) = abs(HTrans(dtftIndexTrans));
    RTrans(notDTFTIndexTrans) = NaN;
    
    % Update the handle structures
    handles.HnumTrans = HnumTrans;
    handles.HdenTrans = HdenTrans;
    handles.RTrans = RTrans; 
    handles.dtftLinIndexTrans = dtftLinIndexTrans;
    handles.notDTFTIndexTrans = notDTFTIndexTrans;
    handles.dtftIndexTrans = dtftIndexTrans;
    handles.HTrans = HTrans;
    handles.w_c = w_c;
    handles.wprime_c = wprime_c;
    handles.wprime_c2 = wprime_c2;
    handles.alpha = alpha;
    handles.beta = beta;
    handles.ZtransNum = ZtransNum; 
    handles.ZtransDen = ZtransDen; 
    
    handles.w = w_1;
    handles.wTransform = wTransform;
    guidata(hObject, handles);
    
end


function [num den pole zero gain] = getInfo(hObject,handles)
% HELP getInfo
%
% Get the correct input and update the other :
% if input = Num/Den => update Zero Pole Gain;
% if input = Zero Pole Gain => update Num/Den; 

if (get(handles.num_den_radiobutton, 'Value'))
    % H = tf(num,den)
    num  = str2num(get(handles.num_edittext, 'String'));
    den  = str2num(get(handles.den_edittext, 'String'));
    % update
    [zero pole gain] = nd2zpg(num, den);
    set(handles.zero_edittext,'String', num2str(zero)); 
    set(handles.pole_edittext,'String', num2str(pole));
    set(handles.gain_edittext,'String', num2str(gain)); 
else
    % H = zpk(zero,pole,1)
    pole = str2num(get(handles.pole_edittext, 'String'));
    zero = str2num(get(handles.zero_edittext, 'String'));
    gain = str2num(get(handles.gain_edittext, 'String'));

    % update
    [num den] = zpg2nd(zero, pole, gain); 
    set(handles.num_edittext,'String', num2str(num)); 
    set(handles.den_edittext,'String', num2str(den)); 
end
guidata(hObject,handles);
end


function [handles] = getH(hObject,handles)
% HELP getH
% Compute the transfer according to the input
% Find the coef of z on the Unit Circle
% 
[Hnum Hden pole zero gain]=getInfo(hObject,handles); 

[x,y] = meshgrid(linspace(-1.25,1.25,1000),linspace(-1.25,1.25,1000)*1j);
z = x+y;
numOrder = length(Hnum) - 1;
denOrder = length(Hden) - 1;

% Calculate H
numSum = 0;
denSum = 0;
    for numPower = 0 : numOrder
        numSum = numSum + Hnum(numPower+1)*z.^(numOrder - numPower); 
    end

    for denPower = 0 : denOrder
        denSum = denSum + Hden(denPower+1)*z.^(denOrder - denPower);
    end

H = numSum./denSum;
% Find coefficients of z that are at the unit circle
dtftLinIndex = find( abs(abs(z) -1) < .001  );
dtftIndex = ind2sub(size(z),dtftLinIndex);
notDTFTIndex = find(abs(z) > 1.025 | abs(z) < .975);
R = zeros(size(z));
R(dtftIndex) = abs(H(dtftIndex));
R(notDTFTIndex) = NaN;


w = angle(z(dtftIndex));
dtft = abs(db(abs(H(dtftIndex))/max(abs(H(dtftIndex)))) +3);
minDTFT = min(abs(db(abs(H(dtftIndex))/max(abs(H(dtftIndex)))) +3)) ;
wcIndex = find(dtft == minDTFT);
wco = w(max(wcIndex));


handles.Hnum = Hnum;
handles.Hden = Hden;
handles.z = z;
handles.x = x;
handles.y = y;
handles.dtftLinIndex = dtftLinIndex;
handles.dtftIndex = dtftIndex;
handles.notDTFTIndex = notDTFTIndex;
handles.R = R; 
handles.H = H;

handles.wco =wco;
guidata(hObject,handles);
end

%%%%%%%%%%%%%Frequency Transform %%%%%% 

% --- Executes on selection change in transform_type_popupmenu.
function transform_type_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to transform_type_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns transform_type_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from transform_type_popupmenu
switch get(handles.transform_type_popupmenu,'Value')
    case 2 % Low Pass
        %fprintf('Im here in low pass\n'); 
       % set(handles.w_c2_edittext,'String', 'NaN'); 
        set(handles.wprime_c2_edittext,'String', 'NaN');
    case 3 % High Pass
        %fprintf('Im here in high pass\n'); 

      %  set(handles.w_c2_edittext,'String', 'NaN'); 
        set(handles.wprime_c2_edittext,'String', 'NaN');
    case 4 % Band Pass
        
    case 5 % Stop Band
        
    otherwise 
end
guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function transform_type_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transform_type_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function w_c_edittext_Callback(hObject, eventdata, handles)
% hObject    handle to w_c_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w_c_edittext as text
%        str2double(get(hObject,'String')) returns contents of w_c_edittext as a double
end

% --- Executes during object creation, after setting all properties.
function w_c_edittext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w_c_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function wprime_c_edittext_Callback(hObject, eventdata, handles)
% hObject    handle to wprime_c_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wprime_c_edittext as text
%        str2double(get(hObject,'String')) returns contents of wprime_c_edittext as a double
end

% --- Executes during object creation, after setting all properties.
function wprime_c_edittext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wprime_c_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function wprime_c2_edittext_Callback(hObject, eventdata, handles)
% hObject    handle to wprime_c2_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wprime_c2_edittext as text
%        str2double(get(hObject,'String')) returns contents of wprime_c2_edittext as a double
end

% --- Executes during object creation, after setting all properties.
function wprime_c2_edittext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wprime_c2_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%%%%PLOT

% --- Executes on button press in Hz_pushbutton.
function Hz_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Hz_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[H x y z dtftLinIndex dtftIndex notDTFTIndex R] = getH(handles);
figure(2);

%figure(1);
%subplot(2,2,1); 

surf(handles.x,real(1j*handles.y),abs(handles.H),'EdgeColor','interp','FaceLighting','none');

    axis([-1.5 1.5 -1.5 1.5 -.5 5])
    caxis([0 7])
    camlight right;
    hold on; 
% Plot the axes

    [x,y] = meshgrid(zeros(1,2), linspace(-1.5,1.5,2) );
    [y1,x1] = meshgrid(zeros(1,2), linspace(-1.5,1.5,2) );
    z = meshgrid(linspace(0,2,2) );
    alphaData = ones(size(z))*0;
    surf(x,y,z,'EdgeColor','none','FaceAlpha',0.5,'FaceColor',[1,.8,.8]);
    surf(x1,y1,z,'EdgeColor','none','FaceAlpha',0.5,'FaceColor',[1,.8,.8]);
    [xCyl,yCyl,zCyl] = cylinder(1,50);
    surf(xCyl,yCyl,zCyl*2,'EdgeColor','interp','EdgeColor','none','FaceAlpha',.5,'FaceColor',[1,.8,.8]);
    hold off;


title('|H(z)|');
xlabel('Real(z)'); 
ylabel('Imag(z)');
zlabel('Magn(H(z))');
end
% --- Executes on button press in Hejw_surf_pushbutton.

function Hejw_surf_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Hejw_surf_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(3); 
%figure(1);
%subplot(2,2,2);
surf(handles.x,real(1j*handles.y),handles.R,'EdgeColor','interp');

title('|H(e\^jw)|');
xlabel('Real(z)'); 
ylabel('Imag(z)');
zlabel('Magn(H(e\^jw))');
end

% --- Executes on button press in Hejw_edittext.0
function Hejw_edittext_Callback(hObject, eventdata, handles)
% hObject    handle to Hejw_edittext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Plot the DTFT (2-dimensions) : |H(e^(jw))| vs w
% For some reason this skips pi * 0.1729 through pi * .8271

figure(4);
%figure(1);
%subplot(2,2,3);
plot(angle(handles.z(handles.dtftIndex))/pi,db(abs(handles.H(handles.dtftIndex))/max(abs(handles.H(handles.dtftIndex)))),'.','MarkerSize',10)
title('|H(e\^jw)|');
xlabel('Frequency [1/pi]'); 
ylabel('Magn(H(e\^jw))');
end

% --- Executes on button press in impulse_pushbutton.
function impulse_pushbutton_Callback(hObject, eventdata, handles)%%

% hObject    handle to impulse_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(5);
%figure(1);
%subplot(2,2,4);
impz(handles.Hnum, handles.Hden); 
end






% --- Executes on button press in Htrans_z_pushbutton.
function Htrans_z_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Htrans_z_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(10);

%figure(1);
%subplot(2,2,1); 

surf(handles.x,real(1j*handles.y),abs(handles.HTrans),'EdgeColor','interp','FaceLighting','none');

    axis([-1.5 1.5 -1.5 1.5 -.5 5])
    caxis([0 7])
    camlight right;
    hold on; 
% Plot the axes

    [x,y] = meshgrid(zeros(1,2), linspace(-1.5,1.5,2) );
    [y1,x1] = meshgrid(zeros(1,2), linspace(-1.5,1.5,2) );
    z = meshgrid(linspace(0,2,2) );
    surf(x,y,z,'EdgeColor','none','FaceAlpha',0.5,'FaceColor',[1,.8,.8]);
    surf(x1,y1,z,'EdgeColor','none','FaceAlpha',0.5,'FaceColor',[1,.8,.8]);
    [xCyl,yCyl,zCyl] = cylinder(1,50);
    surf(xCyl,yCyl,zCyl*2,'EdgeColor','interp','EdgeColor','none','FaceAlpha',.5,'FaceColor',[1,.8,.8]);
    hold off;


title('|HTrans(z)|');
xlabel('Real(z)'); 
ylabel('Imag(z)');
zlabel('Magn(HTrans(z))');

end

% --- Executes on button press in Htrans_ir_pushbutton.
function Htrans_ir_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Htrans_ir_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(9);
%figure(1);
%subplot(2,2,4);
impz(handles.HnumTrans, handles.HdenTrans); 

end

% --- Executes on button press in Htrans_ejw_3D_pushbutton.
function Htrans_ejw_3D_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Htrans_ejw_3D_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(8); 
%figure(1);
%subplot(2,2,2);
surf(handles.x,real(1j*handles.y),handles.RTrans,'EdgeColor','interp');

title('|HTrans(e\^jw)|');
xlabel('Real(z)'); 
ylabel('Imag(z)');
zlabel('Magn(H(e\^jw))');

end

% --- Executes on button press in Htrans_ejw_2D_pushbutton.
function Htrans_ejw_2D_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Htrans_ejw_2D_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(6);
%figure(1);
%subplot(2,2,3);
plot(angle(handles.z(handles.dtftIndex))/pi,db(abs(handles.HTrans(handles.dtftIndex))/max(abs(handles.HTrans(handles.dtftIndex)))),'.','MarkerSize',10)
title('|HTrans(e\^jw)|');
xlabel('Frequency [1/pi]'); 
ylabel('Magn(H(e\^jw))');

end

% --- Executes on button press in allin_H_pushbutton.
function allin_H_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to allin_H_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(16);
title('All In  - H');
end

% --- Executes on button press in allin_Htrans_pushbutton.
function allin_Htrans_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to allin_Htrans_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(17);
title('All In  - HTrans');

end

% --- Executes on button press in freq_trans_vis_pushbutton.
function freq_trans_vis_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to freq_trans_vis_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(20);
title('Frequency transformation Visualisation');
plot(handles.w, handles.wTransform); 

end


