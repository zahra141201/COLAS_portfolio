'use strict';
/** @type {import('sequelize-cli').Migration} */
module.exports = {
  async up(queryInterface, Sequelize) {
    await queryInterface.createTable('Valorations', {
      id: {
        allowNull: false,
        autoIncrement: true,
        primaryKey: true,
        type: Sequelize.INTEGER
      },
      comment: {
        type: Sequelize.STRING
      },
      puntuation: {
        type: Sequelize.INTEGER
      },
      email_user: {
        type: Sequelize.STRING,
        references: { model: 'Users', key: 'email' },
        onUpdate: 'CASCADE',
        onDelete: 'CASCADE'
      },
      made_by: {
        type: Sequelize.STRING,
        references: { model: 'Users', key: 'email' },
        onUpdate: 'CASCADE',
        onDelete: 'CASCADE'
      },
      createdAt: {
        allowNull: false,
        type: Sequelize.DATE
      },
      updatedAt: {
        allowNull: false,
        type: Sequelize.DATE
      }
    });
  },
  async down(queryInterface, Sequelize) {
    await queryInterface.dropTable('Valorations');
  }
};
