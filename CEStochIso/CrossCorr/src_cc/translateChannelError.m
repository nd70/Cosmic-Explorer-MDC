function msg=translateChannelError(chanErr)
%
% translateChannelError --- translate error codes from chanvector and framelist
%                           to the corresponding text description
%
% Chanvector and framelist from the Channel package return a numeric code when
% an error occurs rather than throwing a Matlab error. This routine converts the
% error code to the corresponding text description based on the comments within
% the Channel package.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch (chanErr)
   case 0
    msg = 'No error';
   % Case 1-10 errors from framelist
   case 1
    msg = 'No frame files in list';
   case 2
    msg = 'First frame file missing from disk';
   case 3
    msg = 'Last frame file missing on disk';
   case 4
    msg = 'Interior frame file missing from disk';
   case 5
    msg = 'LIGO_DATAFIND_SERVER not defined';
   % Case 11-30 errors from mllscdatafind
   case 11
    msg = 'gw_data_find calls failed';
   case 12
    msg = 'No frame files returned in list';
   case 13
    msg = 'LIGO_DATAFIND_SERVER not defined';
   % Case 31-40 from fetchseries
   case 31
    msg = 'No files in list';
   case 32
    msg = 'List doesn''t overlap time range';
   case 33
    msg = 'File in list not found';
   case 34
    msg = 'Files do not completely cover time range';
   case 35
    msg = 'Time-series data in files has 0''s';
   case 36
    msg = 'Unable to determine rate';
   case 37
    msg = 'No data for channel';
   % Case 41-50 from chanvector
   case 41
    msg = 'Channel input is not a structure';
   case 42
    msg = 'non-numeric GPS range';
   case 43
    msg = 'Input structure is not a Channel structure';
   case 44
    msg = 'Bad GPS range';
   case 45
    msg = 'External frame list is empty';
   otherwise
    msg = sprintf('Unknown error code %d', chanError);
  end; % switch

return;
