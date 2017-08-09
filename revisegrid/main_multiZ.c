<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='utf-8'>
<title>

GitLab
</title>
<link href="/assets/favicon-6e66d007b09473be55fec600ff017739.ico" rel="shortcut icon" type="image/vnd.microsoft.icon" />
<link href="/assets/application-1bb9d6f4be686c4a4b16debdf17fdf57.css" media="all" rel="stylesheet" />
<link href="/assets/print-8cff922ef0c4bb37621dbac136b43020.css" media="print" rel="stylesheet" />
<script src="/assets/application-f1f8edbae7e79cb56772eed9fa1d4a83.js"></script>
<meta content="authenticity_token" name="csrf-param" />
<meta content="rp7gOGpJ3MqDCwpFgnf0ttSF2HxYvbo4hUIHGrjOiyk=" name="csrf-token" />
<script type="text/javascript">
//<![CDATA[
window.gon={};gon.default_issues_tracker="gitlab";gon.api_version="v3";gon.gravatar_url="http://www.gravatar.com/avatar/%{hash}?s=%{size}\u0026d=mm";gon.relative_url_root="";gon.gravatar_enabled=true;
//]]>
</script>
<meta name="viewport" content="width=device-width, initial-scale=1.0">




</head>

<body class='ui_basic login-page'>
<div class='flash-container'>
<div class='flash-alert'>
You need to sign in before continuing.
</div>
</div>

<div class='container'>
<div class='content'>
<center>
<h1>GitLab</h1>
<p class='light'>
GitLab is open source software to collaborate on code.
<br>
Sign in or browse for <a href="/public/projects">public projects</a>.
</p>
</center>
</div>
</div>
<hr>
<div class='container'>
<div class='content'>
<div class='login-box'>
<h3 class='page-title'>Sign in</h3>
<form accept-charset="UTF-8" action="/users/sign_in" class="new_user" id="new_user" method="post"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /><input name="authenticity_token" type="hidden" value="rp7gOGpJ3MqDCwpFgnf0ttSF2HxYvbo4hUIHGrjOiyk=" /></div><input autofocus="autofocus" class="form-control top" id="user_login" name="user[login]" placeholder="Username or Email" type="text" />
<input class="form-control bottom" id="user_password" name="user[password]" placeholder="Password" type="password" />
<div class='clearfix append-bottom-10'>
<label class='checkbox remember_me' for='user_remember_me'>
<input name="user[remember_me]" type="hidden" value="0" /><input id="user_remember_me" name="user[remember_me]" type="checkbox" value="1" />
<span>Remember me</span>
</label>
</div>
<div>
<input class="btn-create btn" name="commit" type="submit" value="Sign in" />
<div class='pull-right'>
<a class="btn" href="/users/password/new">Forgot your password?</a>
</div>
</div>
</form>



<hr>
<p>
<span class='light'>Did not receive confirmation email?</span>
<a href="/users/confirmation/new">Send again</a>
</p>
</div>

</div>
</div>
</body>
</html>
