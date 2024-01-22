function WarnUser(errorMessage)
	% Alert user via the command window and a popup message.
	fprintf(1, '%s\n', errorMessage); % To command window.
% 	uiwait(warndlg(errorMessage));
	
	% Open the Error Log file for appending.
    cd '/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/'
	fullFileName = 'errr.txt';
	fid = fopen(fullFileName, 'at')
	fprintf(fid, '%s\n', errorMessage); % To file
	fclose(fid);
	return; % from WarnUser()