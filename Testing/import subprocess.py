import subprocess

ffmpeg_path = r"C:\ffmpeg\ffmpeg-2025-06-11-git-f019dd69f0-full_build\bin\ffmpeg.exe"
result = subprocess.run([ffmpeg_path, "-version"], capture_output=True, text=True)
print(result.stdout)
