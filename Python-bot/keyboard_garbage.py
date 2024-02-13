from aiogram.types import ReplyKeyboardRemove, \
    ReplyKeyboardMarkup, KeyboardButton, \
    InlineKeyboardMarkup, InlineKeyboardButton

sort_batteries = 'Сбор_батареек'
sort_plastic_and_glass = "Стекло_и_пластик"
sort_old_clothes = "Старая_одежда"

button_bat = KeyboardButton("/"+sort_batteries) #Текст приветствия, у кнопки обязательный параметр -- текст
button_plas = KeyboardButton("/"+sort_plastic_and_glass)
button_old =  KeyboardButton("/"+sort_old_clothes)

main_kb = ReplyKeyboardMarkup(resize_keyboard=True, one_time_keyboard=True) #разметка клавы для ответа
main_kb.row(button_old, button_plas, button_bat)
#main_kb.add(button_habbits)
#main_kb.add(button_sort)

