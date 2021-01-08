% SENDSMS send message to given phone number(s)
% ** recipient          phone number(s);
% ** subject                subject
% ** msg                message
% Written by YJ based on PSI code 
% *optional*
% ** binary_path        path to sendmail binary file
% ** server             server address
% ** sleep              sleep time between messages
% ** logfile            logfile+path to keep track of past activities
%
% EXAMPLES
%       sendSMS('0041123456789', 'done!');
%       sendSMS({'076123456', '076987654'}, 'done!')
%       sendSMS({'076123456', '076987654'}, {'message1', 'message2'})
%       sendSMS({'076123456', '076987654'}, 'done!', 'sleep', 3600)
%
%

function sendSMS(number, subject, msg)

% =========================================================================
% YOU NEED TO TYPE IN YOUR OWN EMAIL AND PASSWORDS:
mail = 'yijiang.argonne@gmail.com';    %Your GMail email address
password = 'sit_down_reiner';          %Your GMail password
% =========================================================================
%subject = '';
% Format the phone number to 10 digit without dashes
number = strrep(number, '-', '');
if length(number) == 11 && number(1) == '1';
    number = number(2:11);
end
carrier = 'att';
% Information found from
% http://www.sms411.net/2006/07/how-to-send-email-to-phone.html
switch strrep(strrep(lower(carrier),'-',''),'&','')
    case 'alltel';    emailto = strcat(number,'@message.alltel.com');
    case 'att';       emailto = strcat(number,'@txt.att.net');
    case 'boost';     emailto = strcat(number,'@myboostmobile.com');
    case 'cingular';  emailto = strcat(number,'@cingularme.com');
    case 'cingular2'; emailto = strcat(number,'@mobile.mycingular.com');
    case 'nextel';    emailto = strcat(number,'@messaging.nextel.com');
    case 'sprint';    emailto = strcat(number,'@messaging.sprintpcs.com');
    case 'tmobile';   emailto = strcat(number,'@tmomail.net');
    case 'verizon';   emailto = strcat(number,'@vtext.com');
    case 'virgin';    emailto = strcat(number,'@vmobl.com');
end

%% Set up Gmail SMTP service.
% Note: following code found from
% http://www.mathworks.com/support/solutions/data/1-3PRRDV.html
% If you have your own SMTP server, replace it with yours.

% Then this code will set up the preferences properly:
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);

% The following four lines are necessary only if you are using GMail as
% your SMTP server. Delete these lines wif you are using your own SMTP
% server.
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%% Send the email
sendmail(emailto,subject,msg)

end

