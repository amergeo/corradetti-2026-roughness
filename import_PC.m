function [Raw_surface, filename] = import_PC()
% Imports a point cloud (.csv, .xyz, .txt, .ply) and returns it as a
% table with columns {x_col, y_col, z_col}.

[filename, pathname] = uigetfile({'*.csv;*.xyz;*.ply;*.laz;*.txt', ...
    'Point Cloud Files (*.csv, *.xyz, *.ply, *.laz, *.txt)'}, 'Select a Point Cloud File');
fullfilename = strcat(pathname, filename);
[~, ~, ext] = fileparts(filename);

switch ext
    case {'.csv', '.xyz', '.txt'}
        % Preview first 5 lines
        fileID = fopen(fullfilename, 'r');
        previewData = textscan(fileID, '%s', 5, 'Delimiter', '\n', 'Whitespace', '');
        fclose(fileID);
        fprintf('File Preview (first 5 lines):\n');
        for i = 1:min(5, length(previewData{1}))
            disp(previewData{1}{i});
        end

        % Ask for delimiter
        delimiterInput = inputdlg({'Enter the delimiter (space, comma, tab, etc.):'}, ...
            'Delimiter Selection', [1 50], {'tab'});
        delimiter = delimiterInput{1};
        switch lower(delimiter)
            case 'space'; delimiter = ' ';
            case 'comma'; delimiter = ',';
            case 'tab';   delimiter = '\t';
        end

        % Reopen and preview with chosen delimiter
        fileID = fopen(fullfilename, 'r');
        formatSpec = repmat('%s', 1, 10);
        previewData = textscan(fileID, formatSpec, 5, 'Delimiter', delimiter, 'ReturnOnError', false);
        fclose(fileID);
        previewTable = table();
        max_cols = length(previewData);
        for col = 1:max_cols
            previewTable.(sprintf('Col%d', col)) = previewData{col};
        end
        disp('Preview of the first 5 rows (after delimiter selection):');
        disp(previewTable);

        % Ask for column assignment
        colAnswer = inputdlg({'Enter the index of the x column:', ...
            'Enter the index of the y column:', 'Enter the index of the z column:'}, ...
            'Column Assignment', [1 50], {'1', '2', '3'});
        x_col = str2double(colAnswer{1});
        y_col = str2double(colAnswer{2});
        z_col = str2double(colAnswer{3});

        % Read full data
        fileID = fopen(fullfilename, 'r');
        formatSpec = repmat('%f', 1, max_cols);
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
        fclose(fileID);

        Raw_surface = table(dataArray{x_col}, dataArray{y_col}, dataArray{z_col}, ...
            'VariableNames', {'x_col', 'y_col', 'z_col'});

    case '.ply'
        ptCloud = pcread(fullfilename);
        x_coords = double(ptCloud.Location(:,1));
        y_coords = double(ptCloud.Location(:,2));
        z_coords = double(ptCloud.Location(:,3));
        Raw_surface = table(x_coords, y_coords, z_coords, 'VariableNames', {'x_col', 'y_col', 'z_col'});

    case '.laz'
        error('LAZ files are not natively supported by MATLAB. Convert to .las or another format first.');

    otherwise
        error('Unsupported file extension.');
end

disp(head(Raw_surface));
end
