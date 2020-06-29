function sEMG = load_data()
    [filename, pathname] = uigetfile('','Load data');
    data = load(strcat(pathname,filename));
    sEMG = data.signals{1,1};
    sEMG.data = data.signals{1,1}.data;
    sEMG.dt = data.signals{1,1}.dt;
    sEMG.t0 = data.signals{1,1}.t0;
    sEMG.chn_num = data.signals{1, 1}.chn_num;
    %assignin('base','sEMG',sEMG);
end