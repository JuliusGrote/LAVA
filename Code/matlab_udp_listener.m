% MATLAB UDP listener (save as matlab_udp_listener.m)
% Usage: matlab_udp_listener(5555)
function matlab_udp_listener(port)
% port - UDP port to listen on (default 5555)

if nargin < 1 || isempty(port)
    port = 5555;
end

fprintf('Starting MATLAB UDP listener on port %d\n', port);

u = udpport("datagram","LocalPort",port);
cleanupObj = onCleanup(@() clear("u")); 
fprintf('UDP listener started. Waiting for triggers...\n');
triggerReceived = false;
while ~triggerReceived
    if getPendingDatagramCount(u) == 0
        pause(0.0001);
        continue;
    end

    msg = readNextDatagram(u);
    if isempty(msg)
        continue;
    end
    fprintf('[%s] Trigger received: %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'), msg);
    triggerReceived = true;
end
end

function count = getPendingDatagramCount(u)
% getPendingDatagramCount Helper that safely checks for queued data across MATLAB releases
    count = 0;

    if isprop(u,'NumDatagramsAvailable')
        try
            count = double(u.NumDatagramsAvailable);
            return;
        catch %#ok<CTCH>
        end
    end

    fallbackProps = {'NumBytesAvailable','BytesAvailable'};
    for k = 1:numel(fallbackProps)
        prop = fallbackProps{k};
        if isprop(u, prop)
            try
                val = double(u.(prop)); 
                if val > 0
                    count = 1; % treat as at least one datagram pending
                    return;
                end
            catch %#ok<CTCH>
            end
        end
    end
end

function msg = readNextDatagram(u)
% readNextDatagram Read a single datagram and convert it to text safely
    msg = '';
    try
        payload = read(u, 1, "uint8");

        if isa(payload, 'udpport.datagram.Datagram')
            data = payload.Data;
        else
            data = payload;
        end

        if isempty(data)
            return;
        end

        msg = native2unicode(uint8(data(:)'), 'UTF-8');
    catch ME
        warning(ME.identifier, 'Failed to read UDP datagram: %s', ME.message);
    end
end